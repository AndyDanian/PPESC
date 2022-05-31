from libh import *

############# Calculate the electron repulsion integrals ########################
def e2pot(coord, exp, center, lx, ly, lz, output, dalton_normalization):
    """
    Electron repulsion integrals

    Args:
        coord (list): list 2d with coordinates of the atoms
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
        output (int): Output level for integral calculation
        dalton_normalization (bool): it is used the dalton normalization formule

    Return:
        e2pot (array): array 4d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    e2pot: array = np.zeros((total_nprim,total_nprim,total_nprim,total_nprim),dtype=float)

    for i in range(total_nprim):
        for j in range(total_nprim):
            for k in range(total_nprim):
                for l in range(total_nprim):

                    ee: float = electron_repulsion(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        lx[k],
                        ly[k],
                        lz[k],
                        lx[l],
                        ly[l],
                        lz[l],
                        exp[k],
                        exp[l],
                        coord[center[k]][0],
                        coord[center[k]][1],
                        coord[center[k]][2],
                        coord[center[l]][0],
                        coord[center[l]][1],
                        coord[center[l]][2],
                    )

                    e2pot[i,j,k,l] = (
                          normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                        * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                        * normalization(lx[k], ly[k], lz[k], exp[k], dalton_normalization)
                        * normalization(lx[l], ly[l], lz[l], exp[l], dalton_normalization)
                        * ee
                    )

                    if i == j == k == l == 0:
                        print("ee ",ee)
                        print("normalization ",normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization))

    if output > 0:
        print(
            f"\n ***Electron repulsion atomic integrals, time [s]: {time() - start:.6f}"
        )

    return e2pot

if __name__ == "__main__":
    e2_sto2g = e2pot(
        coord = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]],
        exp = [1.3097564, 0.2331360, 1.3097564, 0.2331360],
        center = [0, 0, 1, 1],
        lx = [0, 0, 0, 0],
        ly = [0, 0, 0, 0],
        lz = [0, 0, 0, 0],
        output=11,
        dalton_normalization = False)
    e2_sto3g = e2pot(
        coord = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]],
        exp = [3.4252509, 0.6239137, 0.1688554, 3.4252509, 0.6239137, 0.1688554],
        center = [0, 0, 0, 1, 1, 1],
        lx = [0, 0, 0, 0, 0, 0],
        ly = [0, 0, 0, 0, 0, 0],
        lz = [0, 0, 0, 0, 0, 0],
        output=11,
        dalton_normalization = False)
      
    #print("e2(STO-2G) : ",e2_sto2g)
    print("e2(STO-2G) : ",np.size(e2_sto2g))