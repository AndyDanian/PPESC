from lib1h import *

############# Calculate the potential one body integrals ########################
def nucpot(charge, atom, coord, exp, center, lx, ly, lz, output, dalton_normalization, driver_time):
    """
    Potential integrals

    Args:
        charge (list): list 1d with the charges
        atom (int): atom index
        coord (list): list 2d with coordinates of the atoms
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
        output (int): Output level for integral calculation
        dalton_normalization (bool): it is used the dalton normalization formule
        drive_time (drv_object): Object to manage the time

    Return:
        pot (array): array 1d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    pot: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    count: int = 0
    for i in range(total_nprim):

        for j in range(i, total_nprim):

            veN: float = nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j],
                ly[j],
                lz[j],
                0,
                0,
                0,
                exp[i],
                exp[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                coord[atom][0],
                coord[atom][1],
                coord[atom][2],
            )

            pot[count] = (
                charge[atom]
                * -normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
                * veN
            )

            count += 1
    if output > 10:
        driver_time.add_name_delta_time(
            name = f"Potential Nucleu Atomic Integrals for {atom + 1}-th Atom",
            delta_time = (time() - start)
                )

    return pot
