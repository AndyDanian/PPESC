from lib1h import *


def pangmomp(
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
    Angular momentum between the gradient of the gaussuan functions

            <\nabla \phi| r_g x p | \nabla \phi>     r_g = r - R_g

    This is an anti--symmetric integral

    Args:
        coord (list): list 2d with coordinates of the atoms
        gauge (list): list 1d with gauge coordinates
        magnetic_component (int): magnetic component
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
        pnstcgop (array): array 1d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    pangmomp: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    count: int = 0

    # angmom
    minus_one_column_row: float = 1.0
    if magnetic_component == 0:
        """X Component"""
        mlb: int = 1
        mlc: int = 2
    elif magnetic_component == 1:
        """Y Componente"""
        mlb: int = 0
        mlc: int = 2
        minus_one_column_row: float = -1.0
    elif magnetic_component == 2:
        """Z Componente"""
        mlb: int = 0
        mlc: int = 1
    else:
        raise ValueError(f"***Error\n\n Component not exist: {magnetic_component}")

    dxilidxi: float = 0.0
    dxiljdxi_x: float = 0.0
    ml: list = [lx, ly, lz]
    for i in range(total_nprim):

        for j in range(i, total_nprim):

            # Real{dx_i phi (x_jdx_k - x_kdx_j)_i dx_i phi} =
            # [(ygSkl^0)(dzSmn^0) - (zgSmn^0)(dySkl^0)]dx_idx_jSij^0 =
            # [Skl^1Dmn^1 - Smn^1Dkl^1]Dij^1,1 =
            # [(E_1^kl + YpgE0^kl)(2bE_0^mn+1-l2bE0^mn-1) -
            #  (E_1^mn + ZpgE0^mn)(2bE_0^kl+1-l2bE0^kl-1)]
            # [4abE_0^i+1j+1 - 2ajE_0^i+1j-1 - 2biE_0^i-1j+1 + ijE_0^i-1j-1]
            dxilidxi_x = (
                ml[magnetic_component][i]
                * ml[magnetic_component][j]
                * hermite_coefficient(
                    ml[magnetic_component][i] - 1,
                    ml[magnetic_component][j] - 1,
                    0,
                    coord[center[i]][magnetic_component]
                    - coord[center[j]][magnetic_component],
                    exp[i],
                    exp[j],
                )
                - 2.0
                * exp[i]
                * ml[magnetic_component][j]
                * hermite_coefficient(
                    ml[magnetic_component][i] + 1,
                    ml[magnetic_component][j] - 1,
                    0,
                    coord[center[i]][magnetic_component]
                    - coord[center[j]][magnetic_component],
                    exp[i],
                    exp[j],
                )
                - 2.0
                * exp[j]
                * ml[magnetic_component][i]
                * hermite_coefficient(
                    ml[magnetic_component][i] - 1,
                    ml[magnetic_component][j] + 1,
                    0,
                    coord[center[i]][magnetic_component]
                    - coord[center[j]][magnetic_component],
                    exp[i],
                    exp[j],
                )
                + 4.0
                * exp[j]
                * exp[i]
                * hermite_coefficient(
                    ml[magnetic_component][i] + 1,
                    ml[magnetic_component][j] + 1,
                    0,
                    coord[center[i]][magnetic_component]
                    - coord[center[j]][magnetic_component],
                    exp[i],
                    exp[j],
                )
            )

            # Real{dx_i phi (x_jdx_k - x_kdx_j)_l dx_i phi} =
            # [(ygSkl^0)(dzSmn^0) - (zgSmn^0)(dySkl^0)]Sij^0 =
            dxiljdxi_x = hermite_coefficient(
                ml[magnetic_component][i],
                ml[magnetic_component][j],
                0,
                coord[center[i]][magnetic_component]
                - coord[center[j]][magnetic_component],
                exp[i],
                exp[j],
            )

            dxilidxi: float = 0.0
            dxiljdxi_a: float = 0.0
            dxiljdxi_b: float = 0.0
            lg_sign_firts_second: int = 1.0
            for mli, mlj in [(mlb, mlc), (mlc, mlb)]:
                # Gaussian Center
                Center_Gaussian_xyz: float = (
                    exp[i] * coord[center[i]][mli] + exp[j] * coord[center[j]][mli]
                )
                Center_Gaussian_xyz /= exp[i] + exp[j]

                lg_derivative_term: float = -2.0 * exp[j] * hermite_coefficient(
                    ml[mlj][i],
                    ml[mlj][j] + 1,
                    0,
                    coord[center[i]][mlj] - coord[center[j]][mlj],
                    exp[i],
                    exp[j],
                ) + ml[mlj][j] * hermite_coefficient(
                    ml[mlj][i],
                    ml[mlj][j] - 1,
                    0,
                    coord[center[i]][mlj] - coord[center[j]][mlj],
                    exp[i],
                    exp[j],
                )

                dxilidxi += (
                    lg_sign_firts_second
                    * dxilidxi_x
                    * (
                        (
                            hermite_coefficient(
                                ml[mli][i],
                                ml[mli][j],
                                1,
                                coord[center[i]][mli] - coord[center[j]][mli],
                                exp[i],
                                exp[j],
                            )
                            + (Center_Gaussian_xyz - gauge[mli])
                            * hermite_coefficient(
                                ml[mli][i],
                                ml[mli][j],
                                0,
                                coord[center[i]][mli] - coord[center[j]][mli],
                                exp[i],
                                exp[j],
                            )
                        )
                        * lg_derivative_term
                    )
                )

                #################################################################
                # dy \phi (ydz) dy \phi = -dz \phi (zdy) dz \phi
                dxiljdxi_a += (
                    lg_sign_firts_second
                    * dxiljdxi_x
                    * (
                        4.0
                        * exp[i]
                        * exp[j]
                        * hermite_coefficient(
                            ml[mli][i] + 1,
                            ml[mli][j] + 1,
                            1,
                            coord[center[i]][mli] - coord[center[j]][mli],
                            exp[i],
                            exp[j],
                        )
                        + 4.0
                        * exp[i]
                        * exp[j]
                        * (Center_Gaussian_xyz - gauge[mli])
                        * hermite_coefficient(
                            ml[mli][i] + 1,
                            ml[mli][j] + 1,
                            0,
                            coord[center[i]][mli] - coord[center[j]][mli],
                            exp[i],
                            exp[j],
                        )
                        - 2.0
                        * ml[mli][i]
                        * exp[j]
                        * hermite_coefficient(
                            ml[mli][i] - 1,
                            ml[mli][j] + 1,
                            1,
                            coord[center[i]][mli] - coord[center[j]][mli],
                            exp[i],
                            exp[j],
                        )
                        - 2.0
                        * ml[mli][i]
                        * exp[j]
                        * (Center_Gaussian_xyz - gauge[mli])
                        * hermite_coefficient(
                            ml[mli][i] - 1,
                            ml[mli][j] + 1,
                            0,
                            coord[center[i]][mli] - coord[center[j]][mli],
                            exp[i],
                            exp[j],
                        )
                        - 2.0
                        * exp[i]
                        * ml[mli][j]
                        * hermite_coefficient(
                            ml[mli][i] + 1,
                            ml[mli][j] - 1,
                            1,
                            coord[center[i]][mli] - coord[center[j]][mli],
                            exp[i],
                            exp[j],
                        )
                        - 2.0
                        * exp[i]
                        * ml[mli][j]
                        * (Center_Gaussian_xyz - gauge[mli])
                        * hermite_coefficient(
                            ml[mli][i] + 1,
                            ml[mli][j] - 1,
                            0,
                            coord[center[i]][mli] - coord[center[j]][mli],
                            exp[i],
                            exp[j],
                        )
                        + ml[mli][i]
                        * ml[mli][j]
                        * hermite_coefficient(
                            ml[mli][i] - 1,
                            ml[mli][j] - 1,
                            1,
                            coord[center[i]][mli] - coord[center[j]][mli],
                            exp[i],
                            exp[j],
                        )
                        + ml[mli][i]
                        * ml[mli][j]
                        * (Center_Gaussian_xyz - gauge[mli])
                        * hermite_coefficient(
                            ml[mli][i] - 1,
                            ml[mli][j] - 1,
                            0,
                            coord[center[i]][mli] - coord[center[j]][mli],
                            exp[i],
                            exp[j],
                        )
                    )
                    * lg_derivative_term
                )

                ############################################################
                # dz \phi (ydz) dz \phi = -dy \phi (zdy) dy \phi
                dxiljdxi_b += (
                    lg_sign_firts_second
                    * dxiljdxi_x
                    * (
                        (
                            hermite_coefficient(
                                ml[mli][i],
                                ml[mli][j],
                                1,
                                coord[center[i]][mli] - coord[center[j]][mli],
                                exp[i],
                                exp[j],
                            )
                            + (Center_Gaussian_xyz - gauge[mli])
                            * hermite_coefficient(
                                ml[mli][i],
                                ml[mli][j],
                                0,
                                coord[center[i]][mli] - coord[center[j]][mli],
                                exp[i],
                                exp[j],
                            )
                        )
                        * (
                            ml[mlj][i]
                            * ml[mlj][j]
                            * (ml[mlj][j] - 1)
                            * hermite_coefficient(
                                ml[mlj][i] - 1,
                                ml[mlj][j] - 2,
                                0,
                                coord[center[i]][mlj] - coord[center[j]][mlj],
                                exp[i],
                                exp[j],
                            )
                            - ml[mlj][i]
                            * 2.0
                            * exp[j]
                            * (2.0 * ml[mlj][j] + 1)
                            * hermite_coefficient(
                                ml[mlj][i] - 1,
                                ml[mlj][j],
                                0,
                                coord[center[i]][mlj] - coord[center[j]][mlj],
                                exp[i],
                                exp[j],
                            )
                            + ml[mlj][i]
                            * 4.0
                            * exp[j]
                            * exp[j]
                            * hermite_coefficient(
                                ml[mlj][i] - 1,
                                ml[mlj][j] + 2,
                                0,
                                coord[center[i]][mlj] - coord[center[j]][mlj],
                                exp[i],
                                exp[j],
                            )
                            - 2.0
                            * exp[i]
                            * ml[mlj][j]
                            * (ml[mlj][j] - 1)
                            * hermite_coefficient(
                                ml[mlj][i] + 1,
                                ml[mlj][j] - 2,
                                0,
                                coord[center[i]][mlj] - coord[center[j]][mlj],
                                exp[i],
                                exp[j],
                            )
                            + 4.0
                            * exp[i]
                            * exp[j]
                            * (2.0 * ml[mlj][j] + 1)
                            * hermite_coefficient(
                                ml[mlj][i] + 1,
                                ml[mlj][j],
                                0,
                                coord[center[i]][mlj] - coord[center[j]][mlj],
                                exp[i],
                                exp[j],
                            )
                            - 8.0
                            * exp[i]
                            * exp[j]
                            * exp[j]
                            * hermite_coefficient(
                                ml[mlj][i] + 1,
                                ml[mlj][j] + 2,
                                0,
                                coord[center[i]][mlj] - coord[center[j]][mlj],
                                exp[i],
                                exp[j],
                            )
                        )
                    )
                )
                lg_sign_firts_second = -1.0

            pangmomp[count] = (
                minus_one_column_row
                * normalization(  # (-1)^(column+row), times product sign
                    lx[i], ly[i], lz[i], exp[i], dalton_normalization
                )
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * (dxilidxi + dxiljdxi_a + dxiljdxi_b)
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )
            count += 1

    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Angular Momentum Atomic Integrals for {magnetic_component} Magnetic Component between gradient gaussians",
            delta_time=(time() - start),
        )

    return pangmomp


if __name__ == "__main__":
    # 6-311++G**
    s = pangmomp(
        coord=[[0.0, 0.0, 0.0586476414], [0.0, 0.0, 1.4045523587]],
        exp=[
            33.865,
            5.09479,
            1.15879,
            0.32584,
            0.102741,
            0.036,
            0.75,
            0.75,
            0.75,
            33.865,
            5.09479,
            1.15879,
            0.32584,
            0.102741,
            0.036,
            0.75,
            0.75,
            0.75,
        ],
        center=[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        lx=[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        ly=[0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        lz=[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        magnetic_component=2,
        gauge=[0.0, 0.0, 0.0],
        output=11,
        dalton_normalization=False,
    )

    print("pangmomp : ", s)
