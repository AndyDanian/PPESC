from lib1h import *


def pxaby(
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
    Potential vector of external magnetic field cross 
    lineal momentum

            d/dB_y (A_B x nabla)

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
        curllxp (array): array 1d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    pxaby: list = [0 for i in range(int(total_nprim * total_nprim))]

    count: int = 0

    if magnetic_component == 0:
        """X Componente"""
        mla: list = ly
        mlb: list = lx
        mlc: list = lz
        ra: int = 1
        rb: int = 0
        rc: int = 2
    elif magnetic_component == 1:
        """Y Componente"""
        mla: list = lz
        mlb: list = lx
        mlc: list = ly
        ra: int = 2
        rb: int = 0
        rc: int = 1
    elif magnetic_component == 2:
        """Z Componente"""
        mla: list = ly
        mlb: list = lz
        mlc: list = lx
        ra: int = 1
        rb: int = 2
        rc: int = 0
    elif magnetic_component not in [0, 1, 2]:
        raise ValueError(f"***Error\n\n Component not exist: {magnetic_component}")

    for j in range(total_nprim):

        for i in range(total_nprim):

            sij: list = []
            for ri, mli in enumerate([lx, ly, lz]):
                sij.append(hermite_coefficient(
                                mli[i],
                                mli[j],
                                0,
                                coord[center[i]][ri] - coord[center[j]][ri],
                                exp[i],
                                exp[j],
                                )
                            )

            if magnetic_component == 1:
                ridi: list = []
                for ri, mli in zip([ra, rb], [mla,mlb]):
                    ridi.append(-(2.0*exp[j]*hermite_coefficient(
                                mli[i],
                                mli[j] + 2,
                                0,
                                coord[center[i]][ri] - coord[center[j]][ri],
                                exp[i],
                                exp[j],
                            )
                            -mli[j]*hermite_coefficient(
                                mli[i],
                                mli[j],
                                0,
                                coord[center[i]][ri] - coord[center[j]][ri],
                                exp[i],
                                exp[j],
                            )
                            +(coord[center[j]][ri]-gauge[ri])*
                            (
                            2.0*exp[j]*hermite_coefficient(
                                mli[i],
                                mli[j] + 1,
                                0,
                                coord[center[i]][ri] - coord[center[j]][ri],
                                exp[i],
                                exp[j],
                            )
                            -mli[j]*hermite_coefficient(
                                mli[i],
                                mli[j] - 1,
                                0,
                                coord[center[i]][ri] - coord[center[j]][ri],
                                exp[i],
                                exp[j],
                            )   
                            ))
                        )                    

                temp: float = (sij[rc]*
                                    (ridi[0]*sij[rb]+ridi[1]*sij[ra]
                                    +
                                    2.0*sij[ra]*sij[rb]
                                    )
                                )
            else:
                deri: float = (2.0 * exp[j] * 
                                    hermite_coefficient(
                                    mla[i],
                                    mla[j] + 1,
                                    0,
                                    coord[center[i]][ra] - coord[center[j]][ra],
                                    exp[i],
                                    exp[j],
                                    )
                              - mla[j] *
                                    hermite_coefficient(
                                    mla[i],
                                    mla[j] - 1,
                                    0,
                                    coord[center[i]][ra] - coord[center[j]][ra],
                                    exp[i],
                                    exp[j],
                                    )
                            )
                smn: float = (hermite_coefficient(
                                mlb[i],
                                mlb[j] + 1,
                                0,
                                coord[center[i]][rb] - coord[center[j]][rb],
                                exp[i],
                                exp[j],
                                )
                             + (coord[center[j]][rb]-gauge[rb]) * sij[rb]
                             )
                temp: float = sij[rc]*deri*smn

            pxaby[count] = (
                  normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * temp
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )

            count += 1

    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Curl between nabla with potential vector of magnetic field on y {magnetic_component} ",
            delta_time=(time() - start),
        )

    return pxaby

if __name__ == "__main__":
    # STO-2G
    # print("\n LiH \n")
    # s = pxaby(
    #     coord=[[0.0, 0.0, -0.545857052],[0.0, 0.0, 2.309057052]],
    #     gauge=[0.,0.,0.],
    #     exp=[
    #         6.1638450, 1.0971610, 0.2459160, 0.0623710,
    #         0.2459160, 0.2459160, 0.2459160,
    #         0.0623710, 0.0623710, 0.0623710,
    #         1.3097564, 0.2331360
    #     ],
    #     center=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
    #     lx=[0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0],
    #     ly=[0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0],
    #     lz=[0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0],
    #     magnetic_component=0,
    #     output=9,
    #     dalton_normalization=False,
    #     driver_time=None,
    # )

    # print("p X A_B : ", s, "\n", len(s), "\n\n")

    print(" He STO-2G \n")

    s = pxaby(
        coord=[[0.0, 0.0, 0.0]],
        gauge=[0.,0.,0.],
        exp=[
            2.4328790, 0.4330510
        ],
        center=[0, 0],
        lx=[0, 0],
        ly=[0, 0],
        lz=[0, 0],
        magnetic_component=2,
        output=9,
        dalton_normalization=False,
        driver_time=None,
    )

    print("p X A_B : ", s, "\n", len(s), "\n\n")
