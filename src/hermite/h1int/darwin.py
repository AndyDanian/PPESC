from lib1h import *

def darwin(charge, coord, exp, center, lx, ly, lz, output, dalton_normalization):
    """
    Darwin atomic integrals

    Args:
        charge (list): list 1d with the charges
        coord (list): list 2d with coordinates of the atoms
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

    darwin: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    SPEED_LIGHT: float = 137.0359998
    CONST_ALPHA: float = 1.0 / (SPEED_LIGHT ** 2)
    CONST_DARWIN: float = np.pi * CONST_ALPHA / 2.0

    for k in range(len(coord)):
        count: int = 0

        for i in range(total_nprim):

            for j in range(i, total_nprim):

                dw = gaussian_mult(
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
                    coord[k][0],
                    coord[k][1],
                    coord[k][2],
                )

                darwin[count] += (
                    CONST_DARWIN
                    * charge[k]
                    * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                    * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                    * dw
                )
                count += 1
    if output > 0:
        print_time(
            name = f"Darwin Atomic Integrals", delta_time = (time() - start)
        )

    return darwin
