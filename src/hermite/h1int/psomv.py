from lib1h import *


def psomv(
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
    Mass-Velocity energy correction to the paramagnetic spin-orbit atomic integrals

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
        psomv (array): array 1d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    psomv: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]
    # psomv: list = [0 for i in range(int(total_nprim * total_nprim))]
    c: float = 137.0359998
    alpha: float = 1.0 / c
    ctepsomv: float = alpha * alpha / 8.0

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

        for j in range(i, total_nprim):

            #
            idj_jdi_lap_lap: float = 0.0
            lap_lap_idj_jdi: float = 0.0
            for d2x, d2y, d2z in lap:
                if d2x == 2:
                    lap_l_i: float = float(lx[i])
                    lap_l_j: float = float(lx[j])
                elif d2y == 2:
                    lap_l_i: float = float(ly[i])
                    lap_l_j: float = float(ly[j])
                elif d2z == 2:
                    lap_l_i: float = float(lz[i])
                    lap_l_j: float = float(lz[j])

                for k in range(2):  # run over angular momentum
                    if k == 0:  # (a d_b + b d_a)
                        l_x, l_y, l_z = r_x_c, r_y_c, r_z_c  # l_i     r_i
                        r_x, r_y, r_z = r_x_b, r_y_b, r_z_b
                        sign: float = 1.0
                        l_der_i: float = float(der_l_c[i])
                        l_der_j: float = float(der_l_c[j])
                    else:
                        l_x, l_y, l_z = r_x_b, r_y_b, r_z_b
                        r_x, r_y, r_z = r_x_c, r_y_c, r_z_c
                        sign: float = -1.0
                        l_der_i: float = float(der_l_b[i])
                        l_der_j: float = float(der_l_b[j])
                    #         ##########################
                    #         # ! Terms of Lxyz/rk^3dxx
                    idj_jdi_lap_lap += sign * (
                        16.0
                        * exp[j]
                        * exp[j]
                        * exp[j]
                        * exp[j]
                        * (
                            2.0
                            * exp[i]
                            * nuclear_attraction(
                                lx[i] + l_x,
                                ly[i] + l_y,
                                lz[i] + l_z,
                                lx[j] + d2x * 2,
                                ly[j] + d2y * 2,
                                lz[j] + d2z * 2,
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
                            - l_der_i
                            * nuclear_attraction(
                                lx[i] - l_x,
                                ly[i] - l_y,
                                lz[i] - l_z,
                                lx[j] + d2x * 2,
                                ly[j] + d2y * 2,
                                lz[j] + d2z * 2,
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
                        - 16.0
                        * exp[j]
                        * exp[j]
                        * exp[j]
                        * (2.0 * lap_l_j + 3.0)
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
                            - l_der_i
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
                        + 12.0
                        * exp[j]
                        * exp[j]
                        * (2.0 * lap_l_j * lap_l_j + 2.0 * lap_l_j + 1.0)
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
                            - l_der_i
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
                        - 4.0
                        * exp[j]
                        * (
                            2.0 * lap_l_j * lap_l_j * lap_l_j
                            - 3.0 * lap_l_j * lap_l_j
                            + lap_l_j
                        )
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
                            - l_der_i
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
                        + lap_l_j
                        * (lap_l_j - 1.0)
                        * (lap_l_j - 2.0)
                        * (lap_l_j - 3.0)
                        * (
                            2.0
                            * exp[i]
                            * nuclear_attraction(
                                lx[i] + l_x,
                                ly[i] + l_y,
                                lz[i] + l_z,
                                lx[j] - d2x * 2,
                                ly[j] - d2y * 2,
                                lz[j] - d2z * 2,
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
                            - l_der_i
                            * nuclear_attraction(
                                lx[i] - l_x,
                                ly[i] - l_y,
                                lz[i] - l_z,
                                lx[j] - d2x * 2,
                                ly[j] - d2y * 2,
                                lz[j] - d2z * 2,
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

                    # ! LAP LAP PSO
                    lap_lap_idj_jdi -= sign * (
                        16.0
                        * exp[i]
                        * exp[i]
                        * exp[i]
                        * exp[i]
                        * (
                            2.0
                            * exp[j]
                            * nuclear_attraction(
                                lx[i] + d2x * 2,
                                ly[i] + d2y * 2,
                                lz[i] + d2z * 2,
                                lx[j] + l_x,
                                ly[j] + l_y,
                                lz[j] + l_z,
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
                            - l_der_j
                            * nuclear_attraction(
                                lx[i] + d2x * 2,
                                ly[i] + d2y * 2,
                                lz[i] + d2z * 2,
                                lx[j] - l_x,
                                ly[j] - l_y,
                                lz[j] - l_z,
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
                        - 16.0
                        * exp[i]
                        * exp[i]
                        * exp[i]
                        * (2.0 * lap_l_i + 3.0)
                        * (
                            2.0
                            * exp[j]
                            * nuclear_attraction(
                                lx[i] + d2x,
                                ly[i] + d2y,
                                lz[i] + d2z,
                                lx[j] + l_x,
                                ly[j] + l_y,
                                lz[j] + l_z,
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
                            - l_der_j
                            * nuclear_attraction(
                                lx[i] + d2x,
                                ly[i] + d2y,
                                lz[i] + d2z,
                                lx[j] - l_x,
                                ly[j] - l_y,
                                lz[j] - l_z,
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
                        + 12.0
                        * exp[i]
                        * exp[i]
                        * (2.0 * lap_l_i * lap_l_i + 2.0 * lap_l_i + 1.0)
                        * (
                            2.0
                            * exp[j]
                            * nuclear_attraction(
                                lx[i],
                                ly[i],
                                lz[i],
                                lx[j] + l_x,
                                ly[j] + l_y,
                                lz[j] + l_z,
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
                            - l_der_j
                            * nuclear_attraction(
                                lx[i],
                                ly[i],
                                lz[i],
                                lx[j] - l_x,
                                ly[j] - l_y,
                                lz[j] - l_z,
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
                        - 4.0
                        * exp[i]
                        * (
                            2.0 * lap_l_i * lap_l_i * lap_l_i
                            - 3.0 * lap_l_i * lap_l_i
                            + lap_l_i
                        )
                        * (
                            2.0
                            * exp[j]
                            * nuclear_attraction(
                                lx[i] - d2x,
                                ly[i] - d2y,
                                lz[i] - d2z,
                                lx[j] + l_x,
                                ly[j] + l_y,
                                lz[j] + l_z,
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
                            - l_der_j
                            * nuclear_attraction(
                                lx[i] - d2x,
                                ly[i] - d2y,
                                lz[i] - d2z,
                                lx[j] - l_x,
                                ly[j] - l_y,
                                lz[j] - l_z,
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
                        + lap_l_i
                        * (lap_l_i - 1.0)
                        * (lap_l_i - 2.0)
                        * (lap_l_i - 3.0)
                        * (
                            2.0
                            * exp[j]
                            * nuclear_attraction(
                                lx[i] - d2x * 2,
                                ly[i] - d2y * 2,
                                lz[i] - d2z * 2,
                                lx[j] + l_x,
                                ly[j] + l_y,
                                lz[j] + l_z,
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
                            - l_der_j
                            * nuclear_attraction(
                                lx[i] - d2x * 2,
                                ly[i] - d2y * 2,
                                lz[i] - d2z * 2,
                                lx[j] - l_x,
                                ly[j] - l_y,
                                lz[j] - l_z,
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

            #! Terms cross
            for d2x, d2y, d2z in [(2, 2, 0), (2, 0, 2), (0, 2, 2)]:
                if d2x == 2 and d2y == 2:
                    lap_l_a_i: float = float(lx[i])
                    lap_l_b_i: float = float(ly[i])
                    lap_l_a_j: float = float(lx[j])
                    lap_l_b_j: float = float(ly[j])
                    d2xa, d2ya, d2za = 2, 0, 0
                    d2xb, d2yb, d2zb = 2, -2, 0
                    d2xc, d2yc, d2zc = 0, 2, 0
                    d2xd, d2yd, d2zd = 0, -2, 0
                    d2xe, d2ye, d2ze = -2, 2, 0
                    d2xf, d2yf, d2zf = -2, 0, 0
                elif d2x == 2 and d2z == 2:
                    lap_l_a_i: float = float(lx[i])
                    lap_l_b_i: float = float(lz[i])
                    lap_l_a_j: float = float(lx[j])
                    lap_l_b_j: float = float(lz[j])
                    d2xa, d2ya, d2za = 2, 0, 0
                    d2xb, d2yb, d2zb = 2, 0, -2
                    d2xc, d2yc, d2zc = 0, 0, 2
                    d2xd, d2yd, d2zd = 0, 0, -2
                    d2xe, d2ye, d2ze = -2, 0, 2
                    d2xf, d2yf, d2zf = -2, 0, 0
                elif d2y == 2 and d2z == 2:
                    lap_l_a_i: float = float(ly[i])
                    lap_l_b_i: float = float(lz[i])
                    lap_l_a_j: float = float(ly[j])
                    lap_l_b_j: float = float(lz[j])
                    d2xa, d2ya, d2za = 0, 2, 0
                    d2xb, d2yb, d2zb = 0, 2, -2
                    d2xc, d2yc, d2zc = 0, 0, 2
                    d2xd, d2yd, d2zd = 0, 0, -2
                    d2xe, d2ye, d2ze = 0, -2, 2
                    d2xf, d2yf, d2zf = 0, -2, 0
                for k in range(2):  # run over angular momentum
                    if k == 0:  # (a d_b + b d_a)
                        l_x, l_y, l_z = r_x_c, r_y_c, r_z_c  # l_i     r_i
                        r_x, r_y, r_z = r_x_b, r_y_b, r_z_b
                        sign: float = 1.0
                        l_der_i: float = float(der_l_c[i])
                        l_der_j: float = float(der_l_c[j])
                    else:
                        l_x, l_y, l_z = r_x_b, r_y_b, r_z_b
                        r_x, r_y, r_z = r_x_c, r_y_c, r_z_c
                        sign: float = -1.0
                        l_der_i: float = float(der_l_b[i])
                        l_der_j: float = float(der_l_b[j])
                    ##########################
                    # ! Terms of Lxyz/rk^3dxx
                    idj_jdi_lap_lap += (
                        sign
                        * 2.0
                        * (
                            16.0
                            * exp[j]
                            * exp[j]
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
                                - l_der_i
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
                            + 4.0
                            * exp[j]
                            * exp[j]
                            * (2.0 * lap_l_a_j + 1.0)
                            * (2.0 * lap_l_b_j + 1.0)
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
                                - l_der_i
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
                            + lap_l_a_j
                            * (lap_l_a_j - 1.0)
                            * lap_l_b_j
                            * (lap_l_b_j - 1.0)
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
                                - l_der_i
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
                            #
                            - 8.0
                            * exp[j]
                            * exp[j]
                            * exp[j]
                            * (2.0 * lap_l_b_j + 1.0)
                            * (
                                2.0
                                * exp[i]
                                * nuclear_attraction(
                                    lx[i] + l_x,
                                    ly[i] + l_y,
                                    lz[i] + l_z,
                                    lx[j] + d2xa,
                                    ly[j] + d2ya,
                                    lz[j] + d2za,
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
                                - l_der_i
                                * nuclear_attraction(
                                    lx[i] - l_x,
                                    ly[i] - l_y,
                                    lz[i] - l_z,
                                    lx[j] + d2xa,
                                    ly[j] + d2ya,
                                    lz[j] + d2za,
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
                            + 4.0
                            * exp[j]
                            * exp[j]
                            * lap_l_b_j
                            * (lap_l_b_j - 1.0)
                            * (
                                2.0
                                * exp[i]
                                * nuclear_attraction(
                                    lx[i] + l_x,
                                    ly[i] + l_y,
                                    lz[i] + l_z,
                                    lx[j] + d2xb,
                                    ly[j] + d2yb,
                                    lz[j] + d2zb,
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
                                - l_der_i
                                * nuclear_attraction(
                                    lx[i] - l_x,
                                    ly[i] - l_y,
                                    lz[i] - l_z,
                                    lx[j] + d2xb,
                                    ly[j] + d2yb,
                                    lz[j] + d2zb,
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
                            - 8.0
                            * exp[j]
                            * exp[j]
                            * exp[j]
                            * (2.0 * lap_l_a_j + 1.0)
                            * (
                                2.0
                                * exp[i]
                                * nuclear_attraction(
                                    lx[i] + l_x,
                                    ly[i] + l_y,
                                    lz[i] + l_z,
                                    lx[j] + d2xc,
                                    ly[j] + d2yc,
                                    lz[j] + d2zc,
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
                                - l_der_i
                                * nuclear_attraction(
                                    lx[i] - l_x,
                                    ly[i] - l_y,
                                    lz[i] - l_z,
                                    lx[j] + d2xc,
                                    ly[j] + d2yc,
                                    lz[j] + d2zc,
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
                            * (2.0 * lap_l_a_j + 1.0)
                            * lap_l_b_j
                            * (lap_l_b_j - 1.0)
                            * (
                                2.0
                                * exp[i]
                                * nuclear_attraction(
                                    lx[i] + l_x,
                                    ly[i] + l_y,
                                    lz[i] + l_z,
                                    lx[j] + d2xd,
                                    ly[j] + d2yd,
                                    lz[j] + d2zd,
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
                                - l_der_i
                                * nuclear_attraction(
                                    lx[i] - l_x,
                                    ly[i] - l_y,
                                    lz[i] - l_z,
                                    lx[j] + d2xd,
                                    ly[j] + d2yd,
                                    lz[j] + d2zd,
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
                            + 4.0
                            * exp[j]
                            * exp[j]
                            * lap_l_a_j
                            * (lap_l_a_j - 1.0)
                            * (
                                2.0
                                * exp[i]
                                * nuclear_attraction(
                                    lx[i] + l_x,
                                    ly[i] + l_y,
                                    lz[i] + l_z,
                                    lx[j] + d2xe,
                                    ly[j] + d2ye,
                                    lz[j] + d2ze,
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
                                - l_der_i
                                * nuclear_attraction(
                                    lx[i] - l_x,
                                    ly[i] - l_y,
                                    lz[i] - l_z,
                                    lx[j] + d2xe,
                                    ly[j] + d2ye,
                                    lz[j] + d2ze,
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
                            * (2.0 * lap_l_b_j + 1.0)
                            * lap_l_a_j
                            * (lap_l_a_j - 1.0)
                            * (
                                2.0
                                * exp[i]
                                * nuclear_attraction(
                                    lx[i] + l_x,
                                    ly[i] + l_y,
                                    lz[i] + l_z,
                                    lx[j] + d2xf,
                                    ly[j] + d2yf,
                                    lz[j] + d2zf,
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
                                - l_der_i
                                * nuclear_attraction(
                                    lx[i] - l_x,
                                    ly[i] - l_y,
                                    lz[i] - l_z,
                                    lx[j] + d2xf,
                                    ly[j] + d2yf,
                                    lz[j] + d2zf,
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
                    )

                    #! LAPii LAPjj PSO
                    lap_lap_idj_jdi -= (
                        sign
                        * 2.0
                        * (
                            16.0
                            * exp[i]
                            * exp[i]
                            * exp[i]
                            * exp[i]
                            * (
                                2.0
                                * exp[j]
                                * nuclear_attraction(
                                    lx[i] + d2x,
                                    ly[i] + d2y,
                                    lz[i] + d2z,
                                    lx[j] + l_x,
                                    ly[j] + l_y,
                                    lz[j] + l_z,
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
                                - l_der_j
                                * nuclear_attraction(
                                    lx[i] + d2x,
                                    ly[i] + d2y,
                                    lz[i] + d2z,
                                    lx[j] - l_x,
                                    ly[j] - l_y,
                                    lz[j] - l_z,
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
                            + 4.0
                            * exp[i]
                            * exp[i]
                            * (2.0 * lap_l_a_i + 1.0)
                            * (2.0 * lap_l_b_i + 1.0)
                            * (
                                2.0
                                * exp[j]
                                * nuclear_attraction(
                                    lx[i],
                                    ly[i],
                                    lz[i],
                                    lx[j] + l_x,
                                    ly[j] + l_y,
                                    lz[j] + l_z,
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
                                - l_der_j
                                * nuclear_attraction(
                                    lx[i],
                                    ly[i],
                                    lz[i],
                                    lx[j] - l_x,
                                    ly[j] - l_y,
                                    lz[j] - l_z,
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
                            + lap_l_a_i
                            * (lap_l_a_i - 1.0)
                            * lap_l_b_i
                            * (lap_l_b_i - 1.0)
                            * (
                                2.0
                                * exp[j]
                                * nuclear_attraction(
                                    lx[i] - d2x,
                                    ly[i] - d2y,
                                    lz[i] - d2z,
                                    lx[j] + l_x,
                                    ly[j] + l_y,
                                    lz[j] + l_z,
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
                                - l_der_j
                                * nuclear_attraction(
                                    lx[i] - d2x,
                                    ly[i] - d2y,
                                    lz[i] - d2z,
                                    lx[j] - l_x,
                                    ly[j] - l_y,
                                    lz[j] - l_z,
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
                            #
                            - 8.0
                            * exp[i]
                            * exp[i]
                            * exp[i]
                            * (2.0 * lap_l_b_i + 1.0)
                            * (
                                2.0
                                * exp[j]
                                * nuclear_attraction(
                                    lx[i] + d2xa,
                                    ly[i] + d2ya,
                                    lz[i] + d2za,
                                    lx[j] + l_x,
                                    ly[j] + l_y,
                                    lz[j] + l_z,
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
                                - l_der_j
                                * nuclear_attraction(
                                    lx[i] + d2xa,
                                    ly[i] + d2ya,
                                    lz[i] + d2za,
                                    lx[j] - l_x,
                                    ly[j] - l_y,
                                    lz[j] - l_z,
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
                            + 4.0
                            * exp[i]
                            * exp[i]
                            * lap_l_b_i
                            * (lap_l_b_i - 1.0)
                            * (
                                2.0
                                * exp[j]
                                * nuclear_attraction(
                                    lx[i] + d2xb,
                                    ly[i] + d2yb,
                                    lz[i] + d2zb,
                                    lx[j] + l_x,
                                    ly[j] + l_y,
                                    lz[j] + l_z,
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
                                - l_der_j
                                * nuclear_attraction(
                                    lx[i] + d2xb,
                                    ly[i] + d2yb,
                                    lz[i] + d2zb,
                                    lx[j] - l_x,
                                    ly[j] - l_y,
                                    lz[j] - l_z,
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
                            - 8.0
                            * exp[i]
                            * exp[i]
                            * exp[i]
                            * (2.0 * lap_l_a_i + 1.0)
                            * (
                                2.0
                                * exp[j]
                                * nuclear_attraction(
                                    lx[i] + d2xc,
                                    ly[i] + d2yc,
                                    lz[i] + d2zc,
                                    lx[j] + l_x,
                                    ly[j] + l_y,
                                    lz[j] + l_z,
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
                                - l_der_j
                                * nuclear_attraction(
                                    lx[i] + d2xc,
                                    ly[i] + d2yc,
                                    lz[i] + d2zc,
                                    lx[j] - l_x,
                                    ly[j] - l_y,
                                    lz[j] - l_z,
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
                            * exp[i]
                            * (2.0 * lap_l_a_i + 1.0)
                            * lap_l_b_i
                            * (lap_l_b_i - 1.0)
                            * (
                                2.0
                                * exp[j]
                                * nuclear_attraction(
                                    lx[i] + d2xd,
                                    ly[i] + d2yd,
                                    lz[i] + d2zd,
                                    lx[j] + l_x,
                                    ly[j] + l_y,
                                    lz[j] + l_z,
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
                                - l_der_j
                                * nuclear_attraction(
                                    lx[i] + d2xd,
                                    ly[i] + d2yd,
                                    lz[i] + d2zd,
                                    lx[j] - l_x,
                                    ly[j] - l_y,
                                    lz[j] - l_z,
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
                            + 4.0
                            * exp[i]
                            * exp[i]
                            * lap_l_a_i
                            * (lap_l_a_i - 1.0)
                            * (
                                2.0
                                * exp[j]
                                * nuclear_attraction(
                                    lx[i] + d2xe,
                                    ly[i] + d2ye,
                                    lz[i] + d2ze,
                                    lx[j] + l_x,
                                    ly[j] + l_y,
                                    lz[j] + l_z,
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
                                - l_der_j
                                * nuclear_attraction(
                                    lx[i] + d2xe,
                                    ly[i] + d2ye,
                                    lz[i] + d2ze,
                                    lx[j] - l_x,
                                    ly[j] - l_y,
                                    lz[j] - l_z,
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
                            * exp[i]
                            * (2.0 * lap_l_b_i + 1.0)
                            * lap_l_a_i
                            * (lap_l_a_i - 1.0)
                            * (
                                2.0
                                * exp[j]
                                * nuclear_attraction(
                                    lx[i] + d2xf,
                                    ly[i] + d2yf,
                                    lz[i] + d2zf,
                                    lx[j] + l_x,
                                    ly[j] + l_y,
                                    lz[j] + l_z,
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
                                - l_der_j
                                * nuclear_attraction(
                                    lx[i] + d2xf,
                                    ly[i] + d2yf,
                                    lz[i] + d2zf,
                                    lx[j] - l_x,
                                    ly[j] - l_y,
                                    lz[j] - l_z,
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
                    )

            # NOTE: According to DALTON and WOLFRAM ALPHA, this integral is
            #       PSO LAP LAP
            if lx[i] + ly[i] + lz[i] <= 1 and lx[i] + ly[i] + lz[i] <= 1:
                cte_matematicas: float = 0.5
            else:
                cte_matematicas = 1.0
            psomv[count] = (
                -ctepsomv
                * cte_matematicas
                * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * (lap_lap_idj_jdi + idj_jdi_lap_lap)
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
            )
            # if abs(psomv[count]) > 0.001:
            #     # print("alpha ", exp[i], lx[i], ly[i], lz[i])
            #     # print("beta  ", exp[j], lx[j], ly[j], lz[j])
            #     print(
            #         " LAP LAP PSO ",
            #         -1.0  # ctepsomv
            #         * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
            #         * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
            #         * lap_lap_idj_jdi
            #         * 2.0
            #         * np.pi
            #         / (exp[i] + exp[j]),
            #     )
            #     print(
            #         " PSO LAP LAP  ",
            #         -1.0  # ctepsomv
            #         * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
            #         * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
            #         * idj_jdi_lap_lap
            #         * 2.0
            #         * np.pi
            #         / (exp[i] + exp[j]),
            #     )
            #     print("( ", i + 1, j + 1, " )", psomv[count])
            #     print()
            count += 1
    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Mass-Velocity Energy Correction to the Paramagnetic \
            Spin-Orbit Atomic Integrals, for {spatial_sym} Spatial Symmetry",
            delta_time=(time() - start),
        )

    return psomv


if __name__ == "__main__":
    # STO-2G
    lih: bool = False
    if lih:
        print("\n LiH \n")
        s = psomv(
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
        s = psomv(
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

    print("psomv : ", s, "\n", len(s), "\n\n")
