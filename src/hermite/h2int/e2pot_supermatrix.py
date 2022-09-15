from lib2h import *
from struct import pack
from numba import njit

############# Calculate the electron repulsion integrals ########################
# @njit
def get_values_e2pot(n, coord, exp, center, lx, ly, lz, dalton_normalization):


# (pq|rs) = (qp|rs) = (pq|sr) = (qp|sr) =
# (rs|pq) = (sr|pq) = (rs|qp) = (sr|qp)
    i = j = k = l = 0
    count: int = 0

    with open("AOSPMINT", "wb") as aospmint:
        for i in range(n):
            for j in range(i):
                for k in range(i):

                    if k < i: m: int = k
                    else: m: int = j

                    for l in range(m):

                        # (ij|kl) - 1/4[(ik|jl) + (il|jk)]
                        #  12 34         13 24     14 23
                        ijkl: float = electron_repulsion(
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
                        count += 1
                        # k = j => (ij|kl) = (ij|jl)
                        if k == j:
                            ikjl: float =  ijkl
                        else:
                            ikjl: float = (
                            electron_repulsion(
                            lx[i],
                            ly[i],
                            lz[i],
                            lx[k],
                            ly[k],
                            lz[k],
                            exp[i],
                            exp[k],
                            coord[center[i]][0],
                            coord[center[i]][1],
                            coord[center[i]][2],
                            coord[center[k]][0],
                            coord[center[k]][1],
                            coord[center[k]][2],
                            lx[j],
                            ly[j],
                            lz[j],
                            lx[l],
                            ly[l],
                            lz[l],
                            exp[j],
                            exp[l],
                            coord[center[j]][0],
                            coord[center[j]][1],
                            coord[center[j]][2],
                            coord[center[l]][0],
                            coord[center[l]][1],
                            coord[center[l]][2],
                        ))
                        count += 1
                        if l == j:
                            iljk: float =  ijkl
                        elif l == k:
                            iljk: float = ikjl
                        else:
                            iljk: float = (
                            electron_repulsion(
                            lx[i],
                            ly[i],
                            lz[i],
                            lx[l],
                            ly[l],
                            lz[l],
                            exp[i],
                            exp[l],
                            coord[center[i]][0],
                            coord[center[i]][1],
                            coord[center[i]][2],
                            coord[center[l]][0],
                            coord[center[l]][1],
                            coord[center[l]][2],
                            lx[j],
                            ly[j],
                            lz[j],
                            lx[k],
                            ly[k],
                            lz[k],
                            exp[j],
                            exp[k],
                            coord[center[j]][0],
                            coord[center[j]][1],
                            coord[center[j]][2],
                            coord[center[k]][0],
                            coord[center[k]][1],
                            coord[center[k]][2],
                        )
                        )
                        count += 1

                        #e2pot[i][j][k][l] = (
                        temp_e2pot: float = (
                                normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                            * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                            * normalization(lx[k], ly[k], lz[k], exp[k], dalton_normalization)
                            * normalization(lx[l], ly[l], lz[l], exp[l], dalton_normalization)
                            * (ijkl - 0.25*(ikjl + iljk))
                        )
                        aospmint.write(pack("d",temp_e2pot))

    return count

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

    # e2pot: list = [[[[0.0 for x in range(total_nprim)]
    #                 for x in range(total_nprim)]
    #                 for x in range(total_nprim)]
    #                 for x in range(total_nprim)]

    n: int = total_nprim
    #np.array(e2pot), 
    count = get_values_e2pot(n, np.array(coord), np.array(exp),
                                np.array(center), np.array(lx), np.array(ly), np.array(lz),
                                dalton_normalization)

    if output > 10:
        driver_time.add_name_delta_time(
            name = f"Electron Repulsion Atomic Integrals ({count})", delta_time = (time() - start)
        )

    #return e2pot

if __name__ == "__main__":
    wf = wave_function("../../tests/molden_file/H2_ccpvqz.molden")

    driver_time = drv_time()

    e2pot(
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

    # print("H2 : ",np.size(e2int))
    driver_time.printing()

# H2 CCPVTZ #int: 133672 y tiempo: 0:0:10.070