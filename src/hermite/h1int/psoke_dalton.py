from lib1h import *


def psoke(
    coord,
    spatial_sym,
    atom,
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
    Kinetic energy correction to the paramagnetic spin-orbit atomic integrals

    NOTE: This implementation reproduce values of the DALTON

    Agauges:
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
        drive_time (drv_object): Object to manage the time

    Return:
        psoke (array): array 1d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    psoke: list = [0 for i in range(int(total_nprim * total_nprim))]

    count: int = 0

    r_x_b: int = 0
    r_y_b: int = 0
    r_z_b: int = 0
    r_x_c: int = 0
    r_y_c: int = 0
    r_z_c: int = 0

    if spatial_sym == 0:
        """X Component"""
        r_y_b = 1
        r_z_c = 1
        der_l_b: list = ly
        der_l_c: list = lz
    elif spatial_sym == 1:
        """Y Componente"""
        r_z_b = 1
        r_x_c = 1
        der_l_b: list = lz
        der_l_c: list = lx
    elif spatial_sym == 2:
        """Z Componente"""
        r_x_b = 1
        r_y_c = 1
        der_l_b: list = lx
        der_l_c: list = ly
    else:
        raise ValueError(f"***Error\n\n Component not exist: {spatial_sym}")

    lap: list = [(2, 0, 0), (0, 2, 0), (0, 0, 2)]
    for i in range(total_nprim):

        for j in range(total_nprim):

            #
            idj_jdi_lap: float = 0.0
            for d2x, d2y, d2z in lap:
                if d2x == 2:
                    lap_l_j: float = float(lx[j])
                elif d2y == 2:
                    lap_l_j: float = float(ly[j])
                elif d2z == 2:
                    lap_l_j: float = float(lz[j])

                for k in range(2):
                    if k == 0:  # (a d_b + b d_a)
                        l_x, l_y, l_z = r_x_c, r_y_c, r_z_c  # l_i     r_i
                        r_x, r_y, r_z = r_x_b, r_y_b, r_z_b
                        sign: float = 1.0
                        l_der: float = float(der_l_c[i])
                    else:
                        l_x, l_y, l_z = r_x_b, r_y_b, r_z_b
                        r_x, r_y, r_z = r_x_c, r_y_c, r_z_c
                        sign: float = -1.0
                        l_der: float = float(der_l_b[i])
                    ##########################
                    # ! Terms of Lxyz/rk^3dxx
                    idj_jdi_lap += sign * (
                        4.0
                        * exp[j]
                        * exp[j]
                        * (
                            2.0
                            * exp[i]
                            * nuclear_attraction(
                                lx[i] + l_x,
                                ly[i] + l_y,
                                lz[i] + l_z,
                                lx[j] + d2x,
                                ly[j] + d2y,
                                lz[j] + d2z,
                                r_x,
                                r_y,
                                r_z,
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
                            - l_der
                            * nuclear_attraction(
                                lx[i] - l_x,
                                ly[i] - l_y,
                                lz[i] - l_z,
                                lx[j] + d2x,
                                ly[j] + d2y,
                                lz[j] + d2z,
                                r_x,
                                r_y,
                                r_z,
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
                        )
                        - 2.0
                        * exp[j]
                        * (2.0 * lap_l_j + 1.0)
                        * (
                            2.0
                            * exp[i]
                            * nuclear_attraction(
                                lx[i] + l_x,
                                ly[i] + l_y,
                                lz[i] + l_z,
                                lx[j],
                                ly[j],
                                lz[j],
                                r_x,
                                r_y,
                                r_z,
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
                            - l_der
                            * nuclear_attraction(
                                lx[i] - l_x,
                                ly[i] - l_y,
                                lz[i] - l_z,
                                lx[j],
                                ly[j],
                                lz[j],
                                r_x,
                                r_y,
                                r_z,
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
                        )
                        + lap_l_j
                        * (lap_l_j - 1.0)
                        * (
                            2.0
                            * exp[i]
                            * nuclear_attraction(
                                lx[i] + l_x,
                                ly[i] + l_y,
                                lz[i] + l_z,
                                lx[j] - d2x,
                                ly[j] - d2y,
                                lz[j] - d2z,
                                r_x,
                                r_y,
                                r_z,
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
                            - l_der
                            * nuclear_attraction(
                                lx[i] - l_x,
                                ly[i] - l_y,
                                lz[i] - l_z,
                                lx[j] - d2x,
                                ly[j] - d2y,
                                lz[j] - d2z,
                                r_x,
                                r_y,
                                r_z,
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
                        )
                    )

            # NOTE: According to DALTON and WOLFRAM ALPHA, this integral is
            #       0.5 * PSO LAP
            if lx[i] + ly[i] + lz[i] <= 1 and lx[j] + ly[j] + lz[j] <= 1:
                cte_sp_df = 0.5
            elif (
                lx[i] + ly[i] + lz[i] > 1
                and lx[i] + ly[i] + lz[i] == lx[j] + ly[j] + lz[j]
            ):
                cte_sp_df = 2.0
            else:
                cte_sp_df = 0.5

            #
            psoke[count] = (
                cte_sp_df
                * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * (idj_jdi_lap)
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
            )
            if abs(psoke[count]) > 0.001:
                print(
                    "alpha ",
                    exp[i],
                    lx[i],
                    ly[i],
                    lz[i],
                    normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization),
                )
                print(
                    "beta  ",
                    exp[j],
                    lx[j],
                    ly[j],
                    lz[j],
                    normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization),
                )
                print(count, " : ", psoke[count])
                print()
            count += 1
    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Kinetic-Energy Correction to the Paramagnetic \
            Spin-Orbit Atomic Integrals, for {spatial_sym} Spatial Symmetry",
            delta_time=(time() - start),
        )

    return psoke


if __name__ == "__main__":
    # STO-2G
    lih: bool = False
    if lih:
        print("\n LiH \n")
        s = psoke(
            coord=[[0.0, 0.0, -0.545857052], [0.0, 0.0, 2.309057052]],
            atom=0,
            spatial_sym=0,
            exp=[
                6.1638450,
                1.0971610,
                0.2459160,
                0.0623710,
                0.2459160,
                0.2459160,
                0.2459160,
                0.0623710,
                0.0623710,
                0.0623710,
                1.3097564,
                0.2331360,
            ],
            center=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
            lx=[0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0],
            ly=[0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0],
            lz=[0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0],
            output=9,
            dalton_normalization=False,
            driver_time=None,
        )
    else:
        print("\n He \n")
        s = psoke(
            coord=[[0.0, 0.0, 0.0]],
            atom=0,
            spatial_sym=0,
            exp=[
                9623.91395,
                6.25523565,
                6.25523565,
                6.25523565,
                4.32782104,
                4.32782104,
                4.32782104,
                4.32782104,
                4.32782104,
                4.32782104,
                2.68495795,
                2.68495795,
                2.68495795,
            ],
            center=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            lx=[0, 1, 0, 0, 2, 1, 1, 0, 0, 0, 3, 2, 2],
            ly=[0, 0, 1, 0, 0, 1, 0, 2, 1, 0, 0, 1, 0],
            lz=[0, 0, 0, 1, 0, 0, 1, 0, 1, 2, 0, 0, 1],
            output=9,
            dalton_normalization=True,
            driver_time=None,
        )
    print("psoke : ", s, "\n", len(s), "\n\n")
