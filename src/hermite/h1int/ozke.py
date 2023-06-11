from lib1h import *


def ozke(
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
    Calculates the kinetic energy correction to the orbital Zeeman operator

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
        ozke (array): array 1d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    ozke: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]
    # ozke: list = [0 for i in range(int(total_nprim * (total_nprim)))]

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

            # E_0^ij
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

            # E_1^kl
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

            # ! derecha
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

            # ! derecha
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

            # ! Terms of dxxLx
            xxd = (
                4.0
                * exp[i]
                * exp[i]
                * (
                    hermite_coefficient(
                        spatial_l[i] + 2,
                        spatial_l[j],
                        0,
                        coord[center[i]][magnetic_component]
                        - coord[center[j]][magnetic_component],
                        exp[i],
                        exp[j],
                    )
                )
                - 2.0 * exp[i] * (2.0 * spatial_l[i] + 1.0) * e0ij
                + spatial_l[i]
                * (spatial_l[i] - 1.0)
                * (
                    hermite_coefficient(
                        spatial_l[i] - 2,
                        spatial_l[j],
                        0,
                        coord[center[i]][magnetic_component]
                        - coord[center[j]][magnetic_component],
                        exp[i],
                        exp[j],
                    )
                )
            )

            dxxlx = xxd * ((e1kl + left_rg * e0kl) * pz - (e1mn + right_rg * e0mn) * py)

            # ! Terms of dyyLx
            ########### Dipole
            # E_1^kl + Ypc*E_0^kl
            d_skl_1 = 2.0 * exp[i] * (2.0 * left_l[i] + 1) * (e1kl + left_rg * e0kl)
            # E_1^k-2l + Ypc*E_0^k-2l
            d_skm2l_1 = (
                left_l[i]
                * (left_l[i] - 1)
                * (
                    hermite_coefficient(
                        left_l[i] - 2,
                        left_l[j],
                        1,
                        coord[center[i]][left_coord] - coord[center[j]][left_coord],
                        exp[i],
                        exp[j],
                    )
                    + left_rg
                    * hermite_coefficient(
                        left_l[i] - 2,
                        left_l[j],
                        0,
                        coord[center[i]][left_coord] - coord[center[j]][left_coord],
                        exp[i],
                        exp[j],
                    )
                )
            )
            # E_1^k+2l + Ypc*E_0^k+2l
            d_skt2l_1 = (
                4.0
                * exp[i]
                * exp[i]
                * (
                    hermite_coefficient(
                        left_l[i] + 2,
                        left_l[j],
                        1,
                        coord[center[i]][left_coord] - coord[center[j]][left_coord],
                        exp[i],
                        exp[j],
                    )
                    + left_rg
                    * hermite_coefficient(
                        left_l[i] + 2,
                        left_l[j],
                        0,
                        coord[center[i]][left_coord] - coord[center[j]][left_coord],
                        exp[i],
                        exp[j],
                    )
                )
            )
            ####################### Derivatives
            # E_0^kl+1 - l*E_0^kl-1
            d_dkl_1 = (
                2.0
                * exp[i]
                * (2.0 * left_l[i] + 1.0)
                * (
                    2.0
                    * exp[j]
                    * hermite_coefficient(
                        left_l[i],
                        left_l[j] + 1,
                        0,
                        coord[center[i]][left_coord] - coord[center[j]][left_coord],
                        exp[i],
                        exp[j],
                    )
                    - left_l[j]
                    * hermite_coefficient(
                        left_l[i],
                        left_l[j] - 1,
                        0,
                        coord[center[i]][left_coord] - coord[center[j]][left_coord],
                        exp[i],
                        exp[j],
                    )
                )
            )
            # E_0^k-2l+1 - l*E_0^k-2l-1
            d_dkm2l_1 = (
                left_l[i]
                * (left_l[i] - 1.0)
                * (
                    2.0
                    * exp[j]
                    * hermite_coefficient(
                        left_l[i] - 2,
                        left_l[j] + 1,
                        0,
                        coord[center[i]][left_coord] - coord[center[j]][left_coord],
                        exp[i],
                        exp[j],
                    )
                    - left_l[j]
                    * hermite_coefficient(
                        left_l[i] - 2,
                        left_l[j] - 1,
                        0,
                        coord[center[i]][left_coord] - coord[center[j]][left_coord],
                        exp[i],
                        exp[j],
                    )
                )
            )
            # E_0^k+2l+1 - l*E_0^k+2l-1
            d_dkt2l_1 = (
                4.0
                * exp[i]
                * exp[i]
                * (
                    2.0
                    * exp[j]
                    * hermite_coefficient(
                        left_l[i] + 2,
                        left_l[j] + 1,
                        0,
                        coord[center[i]][left_coord] - coord[center[j]][left_coord],
                        exp[i],
                        exp[j],
                    )
                    - left_l[j]
                    * hermite_coefficient(
                        left_l[i] + 2,
                        left_l[j] - 1,
                        0,
                        coord[center[i]][left_coord] - coord[center[j]][left_coord],
                        exp[i],
                        exp[j],
                    )
                )
            )

            dyylx = (
                pz * (d_skm2l_1 - d_skl_1 + d_skt2l_1)
                + (e1mn + right_rg * e0mn) * (d_dkl_1 - d_dkm2l_1 - d_dkt2l_1)
            ) * e0ij

            # ! Terms of dzzLx
            # E_1^mn + Zpc*E_0^mn
            d_smn_1 = 2.0 * exp[i] * (2.0 * right_l[i] + 1) * (e1mn + right_rg * e0mn)
            # E_1^m-2n + Zpc*E_0^m-2n
            d_smm2n_1 = (
                right_l[i]
                * (right_l[i] - 1)
                * (
                    hermite_coefficient(
                        right_l[i] - 2,
                        right_l[j],
                        1,
                        coord[center[i]][right_coord] - coord[center[j]][right_coord],
                        exp[i],
                        exp[j],
                    )
                    + right_rg
                    * hermite_coefficient(
                        right_l[i] - 2,
                        right_l[j],
                        0,
                        coord[center[i]][right_coord] - coord[center[j]][right_coord],
                        exp[i],
                        exp[j],
                    )
                )
            )
            # E_1^m+2n + Zpc*E_0^m+2n
            d_smt2n_1 = (
                4.0
                * exp[i]
                * exp[i]
                * (
                    hermite_coefficient(
                        right_l[i] + 2,
                        right_l[j],
                        1,
                        coord[center[i]][right_coord] - coord[center[j]][right_coord],
                        exp[i],
                        exp[j],
                    )
                    + right_rg
                    * hermite_coefficient(
                        right_l[i] + 2,
                        right_l[j],
                        0,
                        coord[center[i]][right_coord] - coord[center[j]][right_coord],
                        exp[i],
                        exp[j],
                    )
                )
            )
            ####################### Derivatives
            # E_0^kl+1 - l*E_0^kl-1
            d_dmn_1 = (
                2.0
                * exp[i]
                * (2.0 * right_l[i] + 1.0)
                * (
                    2.0
                    * exp[j]
                    * hermite_coefficient(
                        right_l[i],
                        right_l[j] + 1,
                        0,
                        coord[center[i]][right_coord] - coord[center[j]][right_coord],
                        exp[i],
                        exp[j],
                    )
                    - right_l[j]
                    * hermite_coefficient(
                        right_l[i],
                        right_l[j] - 1,
                        0,
                        coord[center[i]][right_coord] - coord[center[j]][right_coord],
                        exp[i],
                        exp[j],
                    )
                )
            )
            # E_0^k-2l+1 - l*E_0^k-2l-1
            d_dmm2n_1 = (
                right_l[i]
                * (right_l[i] - 1.0)
                * (
                    2.0
                    * exp[j]
                    * hermite_coefficient(
                        right_l[i] - 2,
                        right_l[j] + 1,
                        0,
                        coord[center[i]][right_coord] - coord[center[j]][right_coord],
                        exp[i],
                        exp[j],
                    )
                    - right_l[j]
                    * hermite_coefficient(
                        right_l[i] - 2,
                        right_l[j] - 1,
                        0,
                        coord[center[i]][right_coord] - coord[center[j]][right_coord],
                        exp[i],
                        exp[j],
                    )
                )
            )
            # E_0^k+2l+1 - l*E_0^k+2l-1
            d_dmt2n_1 = (
                4.0
                * exp[i]
                * exp[i]
                * (
                    2.0
                    * exp[j]
                    * hermite_coefficient(
                        right_l[i] + 2,
                        right_l[j] + 1,
                        0,
                        coord[center[i]][right_coord] - coord[center[j]][right_coord],
                        exp[i],
                        exp[j],
                    )
                    - right_l[j]
                    * hermite_coefficient(
                        right_l[i] + 2,
                        right_l[j] - 1,
                        0,
                        coord[center[i]][right_coord] - coord[center[j]][right_coord],
                        exp[i],
                        exp[j],
                    )
                )
            )

            dzzlx = (
                py * (d_smn_1 - d_smm2n_1 - d_smt2n_1)
                + (e1kl + left_rg * e0kl) * (d_dmm2n_1 - d_dmn_1 + d_dmt2n_1)
            ) * e0ij

            # 2 * Lx nabla^2 (signo verificado con matemáticas)
            ozke[count] = (
                # * 0.25 #DALTON (Manninen's Constants)
                -2.0
                * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * (dxxlx + dyylx + dzzlx)
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )
            # if abs(ozke[count]) > 0.01:
            #     print("(", i + 1, j + 1, ") : ", ozke[count])
            count += 1

    if output > 10:
        driver_time.add_name_delta_time(
            "Calculates the Kinetic Energy Correction to the Orbital Zeeman Operator, \
        {magnetic_component} Magnetic Component",
            delta_time=(time() - start),
        )

    return ozke


if __name__ == "__main__":
    # STO-2G
    lih: bool = False
    if lih:
        print("\n LiH \n")
        s = ozke(
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
        s = ozke(
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
    print("ozke : ", s, "\n", len(s), "\n\n")
