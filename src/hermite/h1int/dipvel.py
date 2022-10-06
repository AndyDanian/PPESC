from lib1h import *


def dipvel(
    coord,
    magnetic_component,
    exp,
    center,
    lx,
    ly,
    lz,
    output,
    dalton_normalization,
    driver_time,
):
    """
    Dipole velocity

    Args:
        coord (list): list 2d with coordinates of the atoms
        magnetic_component (int): magnetic component
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
        output (int): Output level for integral calculation
        dalton_normalization (bool): it is used the dalton normalization formule
        drive_time (drv_object): Object to manage the time

    Return:
        dipvel (array): array 1d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    dipvel: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    count: int = 0

    if magnetic_component == 0:
        l_dipvel: list = lx
        l_a: list = ly
        l_b: list = lz
        coord_a: int = 1
        coord_b: int = 2
    elif magnetic_component == 1:
        l_dipvel: list = ly
        l_a: list = lx
        l_b: list = lz
        coord_a: int = 0
        coord_b: int = 2
    elif magnetic_component == 2:
        l_dipvel: list = lz
        l_a: list = ly
        l_b: list = lx
        coord_a: int = 1
        coord_b: int = 0

    for i in range(total_nprim):

        for j in range(i, total_nprim):

            s_a = hermite_coefficient(
                l_a[i],
                l_a[j],
                0,
                coord[center[i]][coord_a] - coord[center[j]][coord_a],
                exp[i],
                exp[j],
            )

            s_b = hermite_coefficient(
                l_b[i],
                l_b[j],
                0,
                coord[center[i]][coord_b] - coord[center[j]][coord_b],
                exp[i],
                exp[j],
            )

            # Horizontal Reccurence
            # int phi d/dx phi dt = Dij^1Skl^0Smn^0 = (2*b*Sij+1^0 - j*Sij-1^0)Skl^0Smn^0
            p_dipvel = 2.0 * exp[j] * hermite_coefficient(
                l_dipvel[i],
                l_dipvel[j] + 1,
                0,
                coord[center[i]][magnetic_component]
                - coord[center[j]][magnetic_component],
                exp[i],
                exp[j],
            ) - l_dipvel[j] * hermite_coefficient(
                l_dipvel[i],
                l_dipvel[j] - 1,
                0,
                coord[center[i]][magnetic_component]
                - coord[center[j]][magnetic_component],
                exp[i],
                exp[j],
            )

            dipvel[count] = (
                normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * p_dipvel
                * s_a
                * s_b
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )
            count += 1
    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Dipole Velocity Atomic Integrals for {magnetic_component} Magnetic Component",
            delta_time=(time() - start),
        )

    return dipvel
