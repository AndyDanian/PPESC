from lib1h import *


def ozmv(
    coord,
    gauge,
    magnetic_component,
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
    Calculates the mass velocity energy correction to the orbital Zeeman operator

    Args:
        coord (list): list 2d with coordinates of the atoms
        gauge (list): list 1d with gauge coordinates
        magnetic_component (int): magnetic component
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
        output (int): Output level for integral calculation
        dalton_normalization (bool): it is used the dalton normalization formule
        drive_time (drv_object): Object to manage the time

    Return:
        ozmv (array): array 1d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    ozmv: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]
    # ozmv: list = [0 for i in range(int(total_nprim * total_nprim))]
    c: float = 137.0359998
    alpha: float = 1.0 / c
    cteozmv: float = alpha * alpha / 8.0

    count: int = 0

    """
    Component Selection L = p x r
                        = (zpy-ypz)x + (xpz-zpx)y + (ypx-xpy)z
    where r = r_e - r_gauge
    """

    if magnetic_component == 0:
        """X Component"""
        left_coord: int = 1
        right_coord: int = 2
        spatial_l: list = lx
        left_l: list = ly
        right_l: list = lz
    elif magnetic_component == 1:
        """Y Componente"""
        left_coord: int = 2
        right_coord: int = 0
        spatial_l: list = ly
        left_l: list = lz
        right_l: list = lx
    elif magnetic_component == 2:
        """Z Componente"""
        left_coord: int = 0
        right_coord: int = 1
        spatial_l: list = lz
        left_l: list = lx
        right_l: list = ly
    else:
        raise ValueError(f"***Error\n\n Component not exist: {magnetic_component}")

    for i in range(total_nprim):

        for j in range(i, total_nprim):

            # Gaussian Center
            left_Pxyz: float = (
                exp[i] * coord[center[i]][left_coord]
                + exp[j] * coord[center[j]][left_coord]
            )
            left_Pxyz = left_Pxyz / (exp[i] + exp[j])
            right_Pxyz: float = (
                exp[i] * coord[center[i]][right_coord]
                + exp[j] * coord[center[j]][right_coord]
            )
            right_Pxyz = right_Pxyz / (exp[i] + exp[j])

            left_rg: float = left_Pxyz - gauge[left_coord]
            right_rg: float = right_Pxyz - gauge[right_coord]

            # & Overlaps ##########################################
            # * E_0^ij = int r_{k,a}^i r_{k,b}^j e^{-(alpha+beta)(r_k^2)} dr_i
            e0ij = hermite_coefficient(
                spatial_l[i],
                spatial_l[j],
                0,
                coord[center[i]][magnetic_component]
                - coord[center[j]][magnetic_component],
                exp[i],
                exp[j],
            )

            # E_0^kl
            e0kl = hermite_coefficient(
                left_l[i],
                left_l[j],
                0,
                coord[center[i]][left_coord] - coord[center[j]][left_coord],
                exp[i],
                exp[j],
            )

            # E_0^mn
            e0mn = hermite_coefficient(
                right_l[i],
                right_l[j],
                0,
                coord[center[i]][right_coord] - coord[center[j]][right_coord],
                exp[i],
                exp[j],
            )
            # & End Overlaps #################################################

            # & int r_{i,a}^k [r_i or r_j] r_{i,b}^l e^{-(alpha+beta)r_i^2} ###
            # * E_1^kl (Gaussian relationship with r_i or r_j)
            e1kl = hermite_coefficient(
                left_l[i],
                left_l[j],
                1,
                coord[center[i]][left_coord] - coord[center[j]][left_coord],
                exp[i],
                exp[j],
            )

            # E_1^mn
            e1mn = hermite_coefficient(
                right_l[i],
                right_l[j],
                1,
                coord[center[i]][right_coord] - coord[center[j]][right_coord],
                exp[i],
                exp[j],
            )
            # & End #############################################################

            # & int r_{i,a}^k [di or dj] r_{i,b}^l e^{-(alpha+beta)r_i^2} ########
            py = 2.0 * exp[j] * hermite_coefficient(
                left_l[i],
                left_l[j] + 1,
                0,
                coord[center[i]][left_coord] - coord[center[j]][left_coord],
                exp[i],
                exp[j],
            ) - left_l[j] * hermite_coefficient(
                left_l[i],
                left_l[j] - 1,
                0,
                coord[center[i]][left_coord] - coord[center[j]][left_coord],
                exp[i],
                exp[j],
            )

            pz = 2.0 * exp[j] * hermite_coefficient(
                right_l[i],
                right_l[j] + 1,
                0,
                coord[center[i]][right_coord] - coord[center[j]][right_coord],
                exp[i],
                exp[j],
            ) - right_l[j] * hermite_coefficient(
                right_l[i],
                right_l[j] - 1,
                0,
                coord[center[i]][right_coord] - coord[center[j]][right_coord],
                exp[i],
                exp[j],
            )
            # & End #################################################################

            # & d4k where k is different to i and j in r_idj - r_jdi ################
            # ! Terms of d4k
            d4k = (
                16.0
                * exp[j]
                * exp[j]
                * exp[j]
                * exp[j]
                * hermite_coefficient(
                    spatial_l[i],
                    spatial_l[j] + 4,
                    0,
                    coord[center[i]][magnetic_component]
                    - coord[center[j]][magnetic_component],
                    exp[i],
                    exp[j],
                )
                - 8.0
                * exp[j]
                * exp[j]
                * exp[j]  #!
                * (4.0 * spatial_l[j] + 6.0)
                * hermite_coefficient(
                    spatial_l[i],
                    spatial_l[j] + 2,
                    0,
                    coord[center[i]][magnetic_component]
                    - coord[center[j]][magnetic_component],
                    exp[i],
                    exp[j],
                )
                + 12.0
                * exp[j]
                * exp[j]  #!
                * (2.0 * spatial_l[j] * spatial_l[j] + 2.0 * spatial_l[j] + 1.0)
                * e0ij
                - 2.0
                * exp[j]
                * spatial_l[j]  #!
                * (spatial_l[j] - 1.0)
                * (4.0 * spatial_l[j] - 2.0)
                * hermite_coefficient(
                    spatial_l[i],
                    spatial_l[j] - 2,
                    0,
                    coord[center[i]][magnetic_component]
                    - coord[center[j]][magnetic_component],
                    exp[i],
                    exp[j],
                )
                + spatial_l[j]
                * (spatial_l[j] - 1.0)
                * (spatial_l[j] - 2.0)  #!
                * (spatial_l[j] - 3.0)
                * hermite_coefficient(
                    spatial_l[i],
                    spatial_l[j] - 4,
                    0,
                    coord[center[i]][magnetic_component]
                    - coord[center[j]][magnetic_component],
                    exp[i],
                    exp[j],
                )
            )

            d4k = d4k * ((e1kl + left_rg * e0kl) * pz - (e1mn + right_rg * e0mn) * py)
            # & End d4k ##################################################################

            # & d4i over idj - jdi #######################################################
            d4i: list = []
            d2i: list = []
            d3i: list = []
            for mli, rgc, xyz, dj, ei, ej in zip(
                [left_l, right_l],
                [left_rg, right_rg],
                [left_coord, right_coord],
                [pz, -py],
                [e1kl + left_rg * e0kl, e1mn + right_rg * e0mn],
                [-(e1mn + right_rg * e0mn), e1kl + left_rg * e0kl],
            ):
                #!
                di4: float = (
                    (
                        16.0
                        * exp[j]
                        * exp[j]
                        * exp[j]
                        * exp[j]
                        * (
                            hermite_coefficient(
                                mli[i],
                                mli[j] + 4,
                                1,
                                coord[center[i]][xyz] - coord[center[j]][xyz],
                                exp[i],
                                exp[j],
                            )
                            + rgc
                            * hermite_coefficient(
                                mli[i],
                                mli[j] + 4,
                                0,
                                coord[center[i]][xyz] - coord[center[j]][xyz],
                                exp[i],
                                exp[j],
                            )
                        )
                    )
                    -
                    # E_1^k+2l + Ypc*E_0^k+2l
                    (
                        8.0
                        * exp[j]
                        * exp[j]
                        * exp[j]
                        * (4.0 * mli[j] + 6.0)
                        * (
                            hermite_coefficient(
                                mli[i],
                                mli[j] + 2,
                                1,
                                coord[center[i]][xyz] - coord[center[j]][xyz],
                                exp[i],
                                exp[j],
                            )
                            + rgc
                            * hermite_coefficient(
                                mli[i],
                                mli[j] + 2,
                                0,
                                coord[center[i]][xyz] - coord[center[j]][xyz],
                                exp[i],
                                exp[j],
                            )
                        )
                    )
                    +
                    # E_1^kl + Ypc*E_0^kl
                    (
                        12.0
                        * exp[j]
                        * exp[j]
                        * (2.0 * mli[j] * mli[j] + 2.0 * mli[j] + 1.0)
                        * ei
                    )
                    -
                    # E_1^k-2l + Ypc*E_0^k-2l
                    (
                        2.0
                        * exp[j]
                        * mli[j]
                        * (mli[j] - 1.0)
                        * (4.0 * mli[j] - 2.0)
                        * (
                            hermite_coefficient(
                                mli[i],
                                mli[j] - 2,
                                1,
                                coord[center[i]][xyz] - coord[center[j]][xyz],
                                exp[i],
                                exp[j],
                            )
                            + rgc
                            * hermite_coefficient(
                                mli[i],
                                mli[j] - 2,
                                0,
                                coord[center[i]][xyz] - coord[center[j]][xyz],
                                exp[i],
                                exp[j],
                            )
                        )
                    )
                    +
                    # E_1^k-4l + Ypc*E_0^k-4l
                    (
                        mli[j]
                        * (mli[j] - 1.0)
                        * (mli[j] - 2.0)
                        * (mli[j] - 3.0)
                        * (
                            hermite_coefficient(
                                mli[i],
                                mli[j] - 4,
                                1,
                                coord[center[i]][xyz] - coord[center[j]][xyz],
                                exp[i],
                                exp[j],
                            )
                            + rgc
                            * hermite_coefficient(
                                mli[i],
                                mli[j] - 4,
                                0,
                                coord[center[i]][xyz] - coord[center[j]][xyz],
                                exp[i],
                                exp[j],
                            )
                        )
                    )
                )
                #!
                # E_0^k+5l+1
                di5: float = (
                    (
                        32.0
                        * exp[j]
                        * exp[j]
                        * exp[j]
                        * exp[j]
                        * exp[j]
                        * hermite_coefficient(
                            mli[i],
                            mli[j] + 5,
                            0,
                            coord[center[i]][xyz] - coord[center[j]][xyz],
                            exp[i],
                            exp[j],
                        )
                    )
                    -
                    # E_0^k+3l+1
                    (
                        16.0
                        * exp[j]
                        * exp[j]
                        * exp[j]
                        * exp[j]
                        * (5.0 * mli[j] + 10.0)
                        * hermite_coefficient(
                            mli[i],
                            mli[j] + 3,
                            0,
                            coord[center[i]][xyz] - coord[center[j]][xyz],
                            exp[i],
                            exp[j],
                        )
                    )
                    +
                    # E_0^kl+1 - l*E_0^kl-1
                    (
                        8.0
                        * exp[j]
                        * exp[j]
                        * exp[j]
                        * (
                            3.0 * (2.0 * mli[j] * mli[j] + 2.0 * mli[j] + 1.0)
                            + (4.0 * mli[j] + 6.0) * (mli[j] + 2.0)
                        )
                        * hermite_coefficient(
                            mli[i],
                            mli[j] + 1,
                            0,
                            coord[center[i]][xyz] - coord[center[j]][xyz],
                            exp[i],
                            exp[j],
                        )
                    )
                    -
                    # E_0^k-1l+1
                    (
                        4.0
                        * exp[j]
                        * exp[j]
                        * mli[j]
                        * (
                            (mli[j] - 1.0) * (4.0 * mli[j] - 2.0)
                            + 3.0 * (2.0 * mli[j] * mli[j] + 2.0 * mli[j] + 1.0)
                        )
                        * hermite_coefficient(
                            mli[i],
                            mli[j] - 1,
                            0,
                            coord[center[i]][xyz] - coord[center[j]][xyz],
                            exp[i],
                            exp[j],
                        )
                    )
                    +
                    # E_0^k-3l+1
                    (
                        2.0
                        * exp[j]
                        * mli[j]
                        * (mli[j] - 1.0)
                        * (mli[j] - 2.0)
                        * (5.0 * mli[j] - 5.0)
                        * hermite_coefficient(
                            mli[i],
                            mli[j] - 3,
                            0,
                            coord[center[i]][xyz] - coord[center[j]][xyz],
                            exp[i],
                            exp[j],
                        )
                    )
                    -
                    # E_0^k-5l+1
                    (
                        mli[j]
                        * (mli[j] - 1.0)
                        * (mli[j] - 2.0)
                        * (mli[j] - 3.0)
                        * (mli[j] - 4.0)
                        * hermite_coefficient(
                            mli[i],
                            mli[j] - 5,
                            0,
                            coord[center[i]][xyz] - coord[center[j]][xyz],
                            exp[i],
                            exp[j],
                        )
                    )
                )
                #! total
                d4i.append((dj * (di4) + ej * (di5)) * e0ij)
                # NOTE: Checked with Wolfram, all right
                # & END ##################################################################

                # & d2id2j r_idj - r_jdi #################################################
                # ! Terms of dyyLx
                d2i.append(
                    # E_1^k+2l + Ypc*E_0^k+2l
                    4.0
                    * exp[j]
                    * exp[j]
                    * (
                        hermite_coefficient(
                            mli[i],
                            mli[j] + 2,
                            1,
                            coord[center[i]][xyz] - coord[center[j]][xyz],
                            exp[i],
                            exp[j],
                        )
                        + rgc
                        * hermite_coefficient(
                            mli[i],
                            mli[j] + 2,
                            0,
                            coord[center[i]][xyz] - coord[center[j]][xyz],
                            exp[i],
                            exp[j],
                        )
                    )
                    # E_1^kl + Ypc*E_0^kl
                    - 2.0 * exp[j] * (2.0 * mli[j] + 1) * ei
                    # E_1^k-2l + Ypc*E_0^k-2l
                    + mli[j]
                    * (mli[j] - 1)
                    * (
                        hermite_coefficient(
                            mli[i],
                            mli[j] - 2,
                            1,
                            coord[center[i]][xyz] - coord[center[j]][xyz],
                            exp[i],
                            exp[j],
                        )
                        + rgc
                        * hermite_coefficient(
                            mli[i],
                            mli[j] - 2,
                            0,
                            coord[center[i]][xyz] - coord[center[j]][xyz],
                            exp[i],
                            exp[j],
                        )
                    )
                )
                #######################
                d3i.append(
                    # E_0^kl+3
                    8.0
                    * exp[j]
                    * exp[j]
                    * exp[j]
                    * hermite_coefficient(
                        mli[i],
                        mli[j] + 3,
                        0,
                        coord[center[i]][xyz] - coord[center[j]][xyz],
                        exp[i],
                        exp[j],
                    )
                    # E_0̣kl+1
                    - 12.0
                    * exp[j]
                    * exp[j]
                    * (mli[j] + 1.0)
                    * hermite_coefficient(
                        mli[i],
                        mli[j] + 1,
                        0,
                        coord[center[i]][xyz] - coord[center[j]][xyz],
                        exp[i],
                        exp[j],
                    )
                    # E_0̣kl-1
                    + 6.0
                    * exp[j]
                    * mli[j]
                    * mli[j]
                    * hermite_coefficient(
                        mli[i],
                        mli[j] - 1,
                        0,
                        coord[center[i]][xyz] - coord[center[j]][xyz],
                        exp[i],
                        exp[j],
                    )
                    # E_0̣kl-3
                    - mli[j]
                    * (mli[j] - 1.0)
                    * (mli[j] - 2.0)
                    * hermite_coefficient(
                        mli[i],
                        mli[j] - 3,
                        0,
                        coord[center[i]][xyz] - coord[center[j]][xyz],
                        exp[i],
                        exp[j],
                    )
                )
            d2id2j: float = -2.0 * (-d2i[0] * d3i[1] + d2i[1] * d3i[0]) * e0ij
            # & END ##################################################################
            # & d2kd2_{i,j} (r_idj - r_jdi) ##########################################
            d2x: float = (
                # E_1^k+2l
                4.0
                * exp[j]
                * exp[j]
                * hermite_coefficient(
                    spatial_l[i],
                    spatial_l[j] + 2,
                    0,
                    coord[center[i]][magnetic_component]
                    - coord[center[j]][magnetic_component],
                    exp[i],
                    exp[j],
                )
                # E_1^kl
                - 2.0 * exp[j] * (2.0 * spatial_l[j] + 1) * e0ij
                # E_1^k-2l
                + spatial_l[j]
                * (spatial_l[j] - 1)
                * hermite_coefficient(
                    spatial_l[i],
                    spatial_l[j] - 2,
                    0,
                    coord[center[i]][magnetic_component]
                    - coord[center[j]][magnetic_component],
                    exp[i],
                    exp[j],
                )
            )
            d2kd2ij: list = []
            d2kd2ij.append(
                2.0 * (pz * d2i[0] - (e1mn + right_rg * e0mn) * d3i[0]) * d2x
            )
            d2kd2ij.append(
                -2.0 * (py * d2i[1] - (e1kl + left_rg * e0kl) * d3i[1]) * d2x
            )
            # & END ##################################################################

            # * L_i nabla^4 nabla^4 ([nabla^4 nabla^4, L_i] = 0)
            ozmv[count] = (
                2.0
                * cteozmv
                * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * (d4k + d4i[0] + d4i[1] + d2kd2ij[0] + d2kd2ij[1] + d2id2j)
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )
            # if abs(ozmv[count]) > 0.01:
            #     print("alpha ", exp[i], lx[i], ly[i], lz[i])
            #     print("beta  ", exp[j], lx[j], ly[j], lz[j])
            #     print(count, " a+1: ", ozmv[count])
            #     print()
            count += 1

    if output > 10:
        driver_time.add_name_delta_time(
            "Calculates the Mass Velocity Energy Correction to the Orbital Zeeman Operator, \
        {magnetic_component} Magnetic Component",
            delta_time=(time() - start),
        )

    return ozmv


if __name__ == "__main__":
    # STO-2G
    lih: bool = False
    if lih:
        print("\n LiH \n")
        s = ozmv(
            coord=[[0.0, 0.0, -0.545857052], [0.0, 0.0, 2.309057052]],
            gauge=[0.0, 0.0, 0.0],
            magnetic_component=0,
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
        s = ozmv(
            coord=[[0.0, 0.0, 0.0]],
            gauge=[0.0, 0.0, 0.0],
            magnetic_component=2,
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
            ],
            center=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            lx=[0, 1, 0, 0, 2, 1, 1, 0, 0, 0],
            ly=[0, 0, 1, 0, 0, 1, 0, 2, 1, 0],
            lz=[0, 0, 0, 1, 0, 0, 1, 0, 1, 2],
            output=9,
            dalton_normalization=True,
            driver_time=None,
        )
    print("ozmv : ", s, "\n", len(s), "\n\n")
