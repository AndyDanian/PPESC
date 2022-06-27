from lib1h import *

def fc(coord, atom, exp, center, lx, ly, lz, output, dalton_normalization):
    """
    Fermi--contact atomic integrals

    Args:
        coord (list): list 2d with coordinates of the atoms
        atom (int): atomic index
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
        output (int): Output level for integral calculation
        dalton_normalization (bool): it is used the dalton normalization formule

    Return:
        angmom (array): array 2d with atomic integrals
    """
    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    fc: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    count: int = 0

    GFACTOR: float = 2.0023193134
    CONST_FC: float = 4.0 * np.pi * GFACTOR / 3.0

    for i in range(total_nprim):

        for j in range(i, total_nprim):

            multiplication_gg: float = gaussian_mult(
                lx[i],
                ly[i],
                lz[i],
                lx[j],
                ly[j],
                lz[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                exp[i],
                exp[j],
                coord[atom][0],
                coord[atom][1],
                coord[atom][2],
            )

            fc[count] = (
                CONST_FC
                * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * multiplication_gg
            )
            count += 1
    if output > 0:
        print_time(
            name = f"Fermi--contact Atomic Integrals for {atom + 1}-th Atom",
            delta_time = (time() - start)
        )

    return fc
