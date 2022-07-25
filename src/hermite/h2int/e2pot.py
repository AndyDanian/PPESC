from lib2h import *
from numba import njit

############# Calculate the electron repulsion integrals ########################
@njit
def get_values_e2pot(e2pot, n, coord, exp, center, lx, ly, lz, dalton_normalization):

    stop: bool = False
    #i = j = k = l = 0
    count: int = 0
    for i in range(n):
        for j in range(i+1):
            for k in range(i+1):
                if k < i: m: int = k + 1
                else: m: int = j + 1
                for l in range(m):

                    calcule: bool = True
                    even_ask_lx: int = (lx[i] + lx[j] + lx[k] + lx[l]) % 2
                    if even_ask_lx != 0 and center[i] == center[j] and center[j] == center[k] and center[k] == center[l]:
                        calcule: bool = False
                    even_ask_ly: int = (ly[i] + ly[j] + ly[k] + ly[l]) % 2
                    if even_ask_ly != 0 and center[i] == center[j] and center[j] == center[k] and center[k] == center[l]:
                        calcule: bool = False
                    even_ask_lz: int = (lz[i] + lz[j] + lz[k] + lz[l]) % 2
                    if even_ask_lz != 0 and center[i] == center[j] and center[j] == center[k] and center[k] == center[l]:
                        calcule: bool = False

                    if abs(e2pot[i][j][k][l]) > 1.0e-12:
                        calcule: bool = False

                    if calcule:
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

                        e2pot[i][j][k][l] = (
                                normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                            * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                            * normalization(lx[k], ly[k], lz[k], exp[k], dalton_normalization)
                            * normalization(lx[l], ly[l], lz[l], exp[l], dalton_normalization)
                            * ee
                        )
                        e2pot[i][j][l][k] = e2pot[j][i][k][l] = e2pot[j][i][l][k] = \
                            e2pot[k][l][i][j] = e2pot[l][k][i][j] = e2pot[k][l][j][i] = e2pot[l][k][j][i] = \
                                e2pot[i][j][k][l]
                        count += 1

    return e2pot, count

def e2pot(coord, exp, center, lx, ly, lz, output, dalton_normalization, driver_time):
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
        drive_time (drv_object): Object to manage the time

    Return:
        e2pot (array): array 4d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    e2pot: list = [[[[0.0 for x in range(total_nprim)]
                    for x in range(total_nprim)]
                    for x in range(total_nprim)]
                    for x in range(total_nprim)]

    n: int = total_nprim #- 1
    e2pot, count = get_values_e2pot(np.array(e2pot), n, np.array(coord), np.array(exp),
                                np.array(center), np.array(lx), np.array(ly), np.array(lz),
                                dalton_normalization)

    if output > 10:
        driver_time.add_name_delta_time(
            name = f"Electron Repulsion Atomic Integrals ({count})", delta_time = (time() - start)
        )

    return e2pot

if __name__ == "__main__":
    wf = wave_function("../../tests/molden_file/H2_ccpvqz.molden")
    ## print wave function informatio
    print("\n Coordinates \n",wf.coordinates)
    print("\n # primtives \n ",wf.primitives_number)
    print("\n Exponents \n",wf.exponents)
    print("\n Centers \n",wf.primitives_centers)
    print("\n mlx \n ",wf.mlx)
    print("\n mly \n ",wf.mly)
    print("\n mlz \n ",wf.mlz)


    driver_time = drv_time()

    e2int = e2pot(
            coord = wf.coordinates,
            exp = wf.exponents,
            center = wf.primitives_centers,
            lx = wf.mlx,
            ly = wf.mly,
            lz = wf.mlz,
            output=11,
            dalton_normalization = False,
            driver_time=driver_time
    )

    print("H2 : ",np.size(e2int))
    driver_time.printing()




    print("\n\n",e2int)
# H2 ccpVTZ Integrals (160831) Time: 0:0:6.737
# H2 ccpVTZ Electron Repulsion Atomic Integrals (136288) Time: 0:0:6.335 con los 4 for