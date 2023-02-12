from lib1h import *


def fcke(coord, atom, exp, center, lx, ly, lz, output, dalton_normalization, driver_time):
    """
    Fermi--contact kinetic atomic integrals

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
        drive_time (drv_object): Object to manage the time

    Return:
        angmom (array): array 2d with atomic integrals
    """
    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    fcke: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]
    count: int = 0

    GFACTOR: float = 2.0023193134
    # CONST_FC: float = 4.0 * np.pi * GFACTOR / 3.0 # DALTON
    CONST_FC: float = 8.0 * np.pi / 3.0

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

            lap_fc: float = 0.0
            fc_lap: float = 0.0
            for lapx, lapy, lapz, mli in zip([2,0,0], [0,2,0], [0,0,2], [lx, ly, lz]):
                lap_fc += (4.0*mli[i]*mli[i]*
                            gaussian_mult(
                                        lx[i] + lapx,
                                        ly[i] + lapy,
                                        lz[i] + lapz,
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
                             -2.0*(2.0*mli[i]+1.0)*multiplication_gg
                            )
                if lx[i] - lapx >= 0 and ly[i] - lapy >= 0 and lz[i] - lapz >= 0:
                    lap_fc +=(mli[i]*(mli[i]-1.0)*
                                gaussian_mult(
                                        lx[i] - lapx,
                                        ly[i] - lapy,
                                        lz[i] - lapz,
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
                            )
                fc_lap += (4.0*mli[j]*mli[j]*
                        gaussian_mult(
                                    lx[i],
                                    ly[i],
                                    lz[i],
                                    lx[j] + lapx,
                                    ly[j] + lapy,
                                    lz[j] + lapz,
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
                        -2.0*(2.0*mli[j]+1.0)*multiplication_gg
                        )
                if lx[j] - lapx >= 0 and ly[j] - lapy >= 0 and lz[j] - lapz >= 0:
                    fc_lap +=(mli[j]*(mli[j]-1.0)*
                                gaussian_mult(
                                    lx[i],
                                    ly[i],
                                    lz[i],
                                    lx[j] - lapx,
                                    ly[j] - lapy,
                                    lz[j] - lapz,
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
                        )

            fcke[count] = (
                CONST_FC
                * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                *(lap_fc  + fc_lap)
            )
            count += 1
    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Fermi--contact Kinetic Atomic Integrals for {atom + 1}-th Atom",
            delta_time=(time() - start),
        )

    return fcke

if __name__ == "__main__":
    # STO-2G
    print("\n LiH \n")
    s = fcke(
        coord=[[0.0, 0.0, -0.545857052],[0.0, 0.0, 2.309057052]],
        atom=0,
        exp=[
            6.1638450, 1.0971610, 0.2459160, 0.0623710,
            0.2459160, 0.2459160, 0.2459160,
            0.0623710, 0.0623710, 0.0623710,
            1.3097564, 0.2331360
        ],
        center=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
        lx=[0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0],
        ly=[0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0],
        lz=[0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0],
        output=9,
        dalton_normalization=False,
        driver_time=None,
    )
    print("p X A_B : ", s, "\n", len(s), "\n\n")