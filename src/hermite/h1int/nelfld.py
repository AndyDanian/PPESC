from lib1h import * 

def nelfld(coord, spatial_sym, atom, exp, center, lx, ly, lz, output, dalton_normalization):
    """
    Spin dipolar atomic integrals, which is a tensor

    Args:
        coord (list): list 2d with coordinates of the atoms
        spatial_sym (int): spatial symmetry index
        atom (int): atomic index
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
        output (int): Output level for integral calculation
        dalton_normalization (bool): it is used the dalton normalization formule

    Return:
        nelfld (array): array 2d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    nelfld: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    count: int = 0

    dx: int = 0
    dy: int = 0
    dz: int = 0
    if spatial_sym == 0:
        dx = 1
    elif spatial_sym == 1:
        dy = 1
    elif spatial_sym == 2:
        dz = 1

    for i in range(total_nprim):

        for j in range(i, total_nprim):

            # Nuclear Electric Field Gradient
            # <phi|xk/rk^3||phi> --> Eqs 9.932, 9.9.18-20

            # *** x
            nef = nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j],
                ly[j],
                lz[j],
                dx,
                dy,
                dz,
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

            nelfld[count] = (
                normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
                * nef
            )
            count += 1

    if output > 0:
        print(
            f"\n *** Nuclear electric field gradient atomic integrals,\
                spatial symmetry {spatial_sym}, time [s]: {time() - start:.6f}"
        )

    return nelfld