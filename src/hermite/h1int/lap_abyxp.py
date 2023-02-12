from lib1h import *


def lap_abyxp(
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

            [Lap, d/dB_y (A_B x nabla)]

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

    lap_abyxp: list = [0 for i in range(int(total_nprim * (total_nprim + 1)/2))]

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
        mla: list = lx
        mlb: list = lz
        mlc: list = ly
        ra: int = 0
        rb: int = 2
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

    for i in range(total_nprim):

        for j in range(i, total_nprim):

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
                lapridi: list = []
                lap: list = []
                for counter, ri, mli in zip([0, 1], [ra, rb], [mla,mlb]):
                    ridi.append(2.0*exp[j]*hermite_coefficient(
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
                            )
                        )
                    #############################################################################
                    temp: float = 0.0
                    for coef, pot in zip([4.0*exp[i]*exp[i],mli[i]*(mli[i] - 1.0)], [2,-2]):
                        temp += coef*(2.0*exp[j]*hermite_coefficient(
                                mli[i] + pot,
                                mli[j] + 2,
                                0,
                                coord[center[i]][ri] - coord[center[j]][ri],
                                exp[i],
                                exp[j],
                            )
                            -mli[j]*hermite_coefficient(
                                mli[i] + pot,
                                mli[j],
                                0,
                                coord[center[i]][ri] - coord[center[j]][ri],
                                exp[i],
                                exp[j],
                            )
                            +(coord[center[j]][ri]-gauge[ri])*
                            (
                            2.0*exp[j]*hermite_coefficient(
                                mli[i] + pot,
                                mli[j] + 1,
                                0,
                                coord[center[i]][ri] - coord[center[j]][ri],
                                exp[i],
                                exp[j],
                            )
                            -mli[j]*hermite_coefficient(
                                mli[i] + pot,
                                mli[j] - 1,
                                0,
                                coord[center[i]][ri] - coord[center[j]][ri],
                                exp[i],
                                exp[j],
                            )   
                            )
                            )
                    lapridi.append(temp-2.0*exp[i]*(2.0*mli[i] + 1.0)*ridi[counter])
                    ###############################################################
                    lap.append(4*exp[i]*exp[i]*
                            hermite_coefficient(
                                mli[i] + 2,
                                mli[j],
                                0,
                                coord[center[i]][ri] - coord[center[j]][ri],
                                exp[i],
                                exp[j],
                            )
                            - 2.0*exp[i]*(2.0*mli[i]+1)*sij[ri]
                            + mli[i]*(mli[i]-1)*
                            hermite_coefficient(
                                mli[i] - 2,
                                mli[j],
                                0,
                                coord[center[i]][ri] - coord[center[j]][ri],
                                exp[i],
                                exp[j],
                            ))
                    counter += 1

                lap.append(4*exp[i]*exp[i]*
                            hermite_coefficient(
                                mlc[i] + 2,
                                mlc[j],
                                0,
                                coord[center[i]][rc] - coord[center[j]][rc],
                                exp[i],
                                exp[j],
                            )
                            - 2.0*exp[i]*(2.0*mlc[i]+1)*sij[rc]
                            + mlc[i]*(mlc[i]-1)*
                            hermite_coefficient(
                                mlc[i] - 2,
                                mlc[j],
                                0,
                                coord[center[i]][rc] - coord[center[j]][rc],
                                exp[i],
                                exp[j],
                            )
                        )

                lapb_abxp: float = (lap[2]*(sij[rb]*ridi[0]+sij[ra]*ridi[1])
                                    +sij[rc]*
                                    (lap[1]*ridi[0]+lapridi[1]*sij[ra]
                                    +lapridi[0]*sij[rb]+lap[0]*ridi[1]
                                    )
                                )

                ### Right term of conmutator [Lap, A_B x P]
                lap: list = []
                for ri, mli in zip([0,1,2], [lx,ly,lz]):
                    lap.append(-(4*exp[j]*exp[j]*
                            hermite_coefficient(
                                mli[i],
                                mli[j] + 2,
                                0,
                                coord[center[i]][ri] - coord[center[j]][ri],
                                exp[i],
                                exp[j],
                            )
                            - 2.0*exp[j]*(2.0*mli[j]+1)*sij[ri]
                            + mli[j]*(mli[j]-1)*
                            hermite_coefficient(
                                mli[i],
                                mli[j] - 2,
                                0,
                                coord[center[i]][ri] - coord[center[j]][ri],
                                exp[i],
                                exp[j],
                            ))
                            )

                deri: list = []
                lapridib: list = []
                for ri, mli in zip([ra,rb],[mla,mlb]):
                    deri.append(
                        2.0*exp[i]*
                        hermite_coefficient(
                                    mli[i] + 2,
                                    mli[j],
                                    0,
                                    coord[center[i]][ri] - coord[center[j]][ri],
                                    exp[i],
                                    exp[j],
                                )
                        -(mli[i]+1.0)*sij[ri]
                        +
                        (coord[center[i]][ri] - gauge[ri])*
                        (
                        2.0*exp[i]*
                        hermite_coefficient(
                                    mli[i] + 1,
                                    mli[j],
                                    0,
                                    coord[center[i]][ri] - coord[center[j]][ri],
                                    exp[i],
                                    exp[j],
                                )
                        - mli[i]*
                        hermite_coefficient(
                                    mli[i] - 1,
                                    mli[j],
                                    0,
                                    coord[center[i]][ri] - coord[center[j]][ri],
                                    exp[i],
                                    exp[j],
                                )
                        )
                    )
                    temp: float = 0.0
                    for poti, coef in zip([2,0,-2], 
                                          [4.0*exp[j]*exp[j], 
                                          -2.0*exp[j]*(2.0*mli[j]+ 1.0),
                                          mli[j]*(mli[j]-1.0)]):
                        temp += (
                            coef*(
                                2.0*exp[i]*
                                hermite_coefficient(
                                            mli[i] + 2,
                                            mli[j] + poti,
                                            0,
                                            coord[center[i]][ri] - coord[center[j]][ri],
                                            exp[i],
                                            exp[j],
                                        )
                                -(mli[i]+1.0)*
                                hermite_coefficient(
                                            mli[i],
                                            mli[j] + poti,
                                            0,
                                            coord[center[i]][ri] - coord[center[j]][ri],
                                            exp[i],
                                            exp[j],
                                        )
                                +
                                (coord[center[i]][ri] - gauge[ri])*
                                (
                                2.0*exp[i]*
                                hermite_coefficient(
                                            mli[i] + 1,
                                            mli[j] + poti,
                                            0,
                                            coord[center[i]][ri] - coord[center[j]][ri],
                                            exp[i],
                                            exp[j],
                                        )
                                - mli[i]*
                                hermite_coefficient(
                                            mli[i] - 1,
                                            mli[j] + poti,
                                            0,
                                            coord[center[i]][ri] - coord[center[j]][ri],
                                            exp[i],
                                            exp[j],
                                        )
                                )
                            )
                        )
                    lapridib.append(temp)

                abxp_lapb: float = (-lap[rc]*(sij[rb]*deri[0]+sij[ra]*deri[1])
                                    +
                                    sij[rc]*(sij[ra]*lapridib[1]-lap[rb]*deri[0])
                                    +
                                    sij[rc]*(lapridib[0]*sij[rb]-lap[ra]*deri[1])
                                    )


                lap_abyxp[count] = (
                    - normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                    * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                    * (lapb_abxp + abxp_lapb)
                    * np.power(np.pi / (exp[i] + exp[j]), 1.5)
                )
            else:
                derij: float = 1.0
                for ri, mli in zip([ra,rb],[mla, mlb]):
                    derij *= (2.0 * exp[i] * hermite_coefficient(
                                mli[i] + 1,
                                mli[j],
                                0,
                                coord[center[i]][ri] - coord[center[j]][ri],
                                exp[i],
                                exp[j],
                                )
                            - mli[i]
                            * hermite_coefficient(
                                mli[i] - 1,
                                mli[j],
                                0,
                                coord[center[i]][ri] - coord[center[j]][ri],
                                exp[i],
                                exp[j],
                                ))

                lap_abyxp[count] = (
                    - normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                    * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                    * 2.0*derij*sij[rc]
                    * np.power(np.pi / (exp[i] + exp[j]), 1.5)
                )
            
            count += 1

    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Curl between potential vector of magnetic field with nabla on {magnetic_component} ",
            delta_time=(time() - start),
        )

    return lap_abyxp

if __name__ == "__main__":
    # STO-2G
    s = lap_abyxp(
        coord=[[0.0, 0.0, -0.545857052],[0.0, 0.0, 2.309057052]],
        gauge=[0.,0.,0.],
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
        magnetic_component=1,
        output=9,
        dalton_normalization=False,
        driver_time=None,
    )

    print("A_B X p : ", s, "\n", len(s))
