from lib1h import *

def diplen(coord, magnetic_component, rdipole, exp, center, lx, ly, lz, output, dalton_normalization, driver_time):
    """
    Dipole lenght atomic integrals

    Args:
        coord (list): list 2d with coordinates of the atoms
        magnetic_component (int): magnetic component
        rdiple (list): list 1d with reference coordinates for the dipole
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
        output (int): Output level for integral calculation
        dalton_normalization (bool): it is used the dalton normalization formule
        drive_time (drv_object): Object to manage the time

    Return:
        diplen (array): array 1d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    diplen: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    count: int = 0

    if magnetic_component == 0:
        l_diplen: list = lx
        l_a: list = ly
        l_b: list = lz
        coord_a: int = 1
        coord_b: int = 2
    elif magnetic_component == 1:
        l_diplen: list = ly
        l_a: list = lx
        l_b: list = lz
        coord_a: int = 0
        coord_b: int = 2
    elif magnetic_component == 2:
        l_diplen: list = lz
        l_a: list = ly
        l_b: list = lx
        coord_a: int = 1
        coord_b: int = 0

    for i in range(total_nprim):

        for j in range(i, total_nprim):

            s_dipole = E(
                l_diplen[i],
                l_diplen[j],
                0,
                coord[center[i]][magnetic_component] - coord[center[j]][magnetic_component], #0
                exp[i],
                exp[j],
            )

            s_a = E(
                l_a[i],
                l_a[j],
                0,
                coord[center[i]][coord_a] - coord[center[j]][coord_a],
                exp[i],
                exp[j],
            )

            s_b = E(
                l_b[i],
                l_b[j],
                0,
                coord[center[i]][coord_b] - coord[center[j]][coord_b],
                exp[i],
                exp[j],
            )
            # Eq 9.5.43 Helgaker
            # <phi|x_c|phi> = <phi|x_p + Xpc|phi> =
            # Sij^1Skl^0Smn^0=(E1^ij + Xpc*E0^ij)E0^klE0^mn(pi/p)^1.5

            Pxyz = (
                exp[i] * coord[center[i]][magnetic_component] #0
                + exp[j] * coord[center[j]][magnetic_component]
            )
            Pxyz = Pxyz / (exp[i] + exp[j])
            rpk = Pxyz - rdipole[magnetic_component] #0

            xyzdipole = E(
                l_diplen[i],
                l_diplen[j],
                1,
                coord[center[i]][magnetic_component] - coord[center[j]][magnetic_component], #0
                exp[i],
                exp[j],
            )

            diplen[count] = (
                normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * (xyzdipole + rpk * s_dipole)
                * s_a
                * s_b
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )
            count += 1
    if output > 0:
        driver_time.add_name_delta_time(
            name = f"Dipole Lenght Atomic Integrals for {magnetic_component} Magnetic Component",
            delta_time = (time() - start)
        )

    return diplen