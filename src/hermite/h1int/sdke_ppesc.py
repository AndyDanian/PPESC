from lib1h import *


def sdke_ppesc(
    coord,
    magnetic_component,
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
    Kinetic energy correction to spin dipolar atomic integrals, which is a tensor

    Args:
        coord (list): list 2d with coordinates of the atoms
        magnetic_component (int): magnetic component
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
        sdke_ppesc (array): array 2d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    # sdke_ppesc: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]
    sdke_ppesc: list = [0 for i in range(int(total_nprim * (total_nprim)))]

    count: int = 0

    GFACTOR: float = 2.0023193134
    # CONST_SDk: float = GFACTOR / 2.0 * 1 / 3.0 # DALTON
    CONST_SDk: float = (
        1 / 3.0
    )  # Due to take 1/r⁵ like derivatives: x_k/r⁵ = 1/3 d²x 1/r

    dx_x: int = 0
    dy_x: int = 0
    dz_x: int = 0
    dx_y: int = 0
    dy_y: int = 0
    dz_y: int = 0
    dx_z: int = 0
    dy_z: int = 0
    dz_z: int = 0

    CTE_XR: float = 1.0
    CTE_YR: float = 1.0
    CTE_ZR: float = 1.0

    x_active: bool = False
    y_active: bool = False
    z_active: bool = False

    if magnetic_component == 0:
        if spatial_sym == 0:
            CTE_XR = 2.0
            CTE_YR = -1.0
            CTE_ZR = -1.0
            dx_x = 2
            dy_y = 2
            dz_z = 2
            x_active = True
            y_active = True
            z_active = True
        elif spatial_sym == 1:
            CTE_YR = 3.0
            y_active = True
            dx_y = 1
            dy_y = 1
        elif spatial_sym == 2:
            CTE_ZR = 3.0
            z_active = True
            dx_z = 1
            dz_z = 1
    if magnetic_component == 1:
        if spatial_sym == 0:
            CTE_XR = 3.0
            x_active = True
            dx_x = 1
            dy_x = 1
        elif spatial_sym == 1:
            CTE_XR = -1.0
            CTE_YR = 2.0
            CTE_ZR = -1.0
            dx_x = 2
            dy_y = 2
            dz_z = 2
            x_active = True
            y_active = True
            z_active = True
        elif spatial_sym == 2:
            CTE_ZR = 3.0
            z_active = True
            dy_z = 1
            dz_z = 1
    if magnetic_component == 2:
        if spatial_sym == 0:
            CTE_XR = 3.0
            x_active = True
            dx_x = 1
            dz_x = 1
        elif spatial_sym == 1:
            CTE_YR = 3.0
            y_active = True
            dy_y = 1
            dz_y = 1
        elif spatial_sym == 2:
            CTE_YR = -1.0
            CTE_XR = -1.0
            CTE_ZR = 2.0
            dx_x = 2
            dy_y = 2
            dz_z = 2
            x_active = True
            y_active = True
            z_active = True

    laplacian: list = [
        (lx, 2, 0, 0),
        (ly, 0, 2, 0),
        (lz, 0, 0, 2),
    ]

    for i in range(total_nprim):

        for j in range(total_nprim):

            # SDk = <phi|(3r_krk^T-r_k^2)/rk^5|phi>
            # https://www.wolframalpha.com/input?i2d=true&i=D%5BDivide%5B1%2CSqrt%5BPower%5Bx%2C2%5D%2BPower%5By%2C2%5D%2BPower%5Bz%2C2%5D%5D%5D%2C%7Bx%2C2%7D%5D
            # !Note: is add the constant 1/3, view the before wolfram result
            lapsd: float = 0.0
            sdlap: float = 0.0
            for rtrue, cter, d2ij, d2kl, d2nm in [
                (x_active, CTE_XR, dx_x, dy_x, dz_x),
                (y_active, CTE_YR, dx_y, dy_y, dz_y),
                (z_active, CTE_ZR, dx_z, dy_z, dz_z),
            ]:
                if rtrue:
                    sdr = nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        d2ij,
                        d2kl,
                        d2nm,
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

                    if lx[i] + ly[i] + lz[i] == 0 or lx[j] + ly[j] + lz[j] == 0:
                        quarter: float = 1.0 / 4.0
                    else:
                        quarter = 1.0
                    for ml, d2x, d2y, d2z in laplacian:

                        lapsd += (
                            quarter
                            * cter
                            * (
                                4.0
                                * exp[i]
                                * exp[i]
                                * nuclear_attraction(
                                    lx[i] + d2x,
                                    ly[i] + d2y,
                                    lz[i] + d2z,
                                    lx[j],
                                    ly[j],
                                    lz[j],
                                    d2ij,
                                    d2kl,
                                    d2nm,
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
                                - 2.0 * exp[i] * (2.0 * ml[i] + 1.0) * sdr
                                + ml[i]
                                * (ml[i] - 1.0)
                                * nuclear_attraction(
                                    lx[i] - d2x,
                                    ly[i] - d2y,
                                    lz[i] - d2z,
                                    lx[j],
                                    ly[j],
                                    lz[j],
                                    d2ij,
                                    d2kl,
                                    d2nm,
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

                        sdlap += (
                            quarter
                            * cter
                            * (
                                4.0
                                * exp[j]
                                * exp[j]
                                * nuclear_attraction(
                                    lx[i],
                                    ly[i],
                                    lz[i],
                                    lx[j] + d2x,
                                    ly[j] + d2y,
                                    lz[j] + d2z,
                                    d2ij,
                                    d2kl,
                                    d2nm,
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
                                - 2.0 * exp[j] * (2.0 * ml[j] + 1.0) * sdr
                                + ml[j]
                                * (ml[j] - 1.0)
                                * nuclear_attraction(
                                    lx[i],
                                    ly[i],
                                    lz[i],
                                    lx[j] - d2x,
                                    ly[j] - d2y,
                                    lz[j] - d2z,
                                    d2ij,
                                    d2kl,
                                    d2nm,
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
            cte_amo_s: float = 1.0
            if magnetic_component != spatial_sym:
                if lx[i] + ly[i] + lz[i] == 0 or lx[j] + ly[j] + lz[j] == 0:
                    cte_amo_s = 1.0 / 4.0

            sdke_ppesc[count] = (
                CONST_SDk
                * cte_amo_s
                * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
                * (lapsd + 3.0 * sdlap)
            )

            # if abs(sdke_ppesc[count]) > 0:
            #     print(
            #         "(",
            #         i + 1,
            #         j + 1,
            #         "): ",
            #         sdke_ppesc[count],
            #     )

            count += 1

    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Spin-Dipolar Atomic Integrals, {spatial_sym} Spatial Symmetry and \
        {magnetic_component} Magnetic Component of {atom + 1}-th Atom",
            delta_time=(time() - start),
        )

    return sdke_ppesc


if __name__ == "__main__":
    # STO-2G
    lih: bool = False
    if lih:
        print("\n LiH \n")
        s = sdke_ppesc(
            coord=[[0.0, 0.0, -0.545857052], [0.0, 0.0, 2.309057052]],
            magnetic_component=0,
            spatial_sym=0,
            atom=0,
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
        s = sdke_ppesc(
            coord=[[0.0, 0.0, 0.0]],
            magnetic_component=2,
            spatial_sym=1,
            atom=0,
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
    print("sdke_ppesc : ", s, "\n", len(s), "\n\n")
