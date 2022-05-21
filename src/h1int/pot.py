from libh import *

############# Calculate the potential one body integrals ########################
def pot(charge, atom, coord, exp, center, lx, ly, lz, output):
    """
    Potential integrals

    Args:
        charge (list): list 1d with the charges
        atom (list): list 1d with atoms index
        coord (list): list 2d with coordinates of the atoms
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
        output (int): Output level for integral calculation

    Return:
        pot (array): array 1d with atomic integrals
    """

    start = time()
    # Primitive total in the cluster
    total_nprim = len(exp)

    pot = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    count = 0
    for i in range(total_nprim):

        for j in range(i, total_nprim):

            veN = nuclear_attraction(
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
                * -Norm[lx[i] + ly[i] + lz[i]](exp[i])
                * Norm[lx[j] + ly[j] + lz[j]](exp[j])
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
                * veN
            )

            count += 1
    if output > 0:
        print(
            f"\n *** One body atomic potential integrals for {atom + 1}° atom, time [s]: {time() - start:.6f}"
        )

    return pot
