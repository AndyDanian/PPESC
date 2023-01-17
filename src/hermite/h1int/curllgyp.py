from lib1h import *


def curllgyp(
    coord,
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
    Curl Angular moment on y by lineal momentum, which is a vector

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

    curllgyp: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    count: int = 0

    """
    Component Selection p x (Lg,y p) = 
                        = dy(zgdx-xgdz)dz-dz(zgdx-xgdz)dy + ...
    where yg = y_e - y_gauge
    """


    if magnetic_component == 0:
        """X Componente"""
        mli: list = ly
        mlj: list = lx
        mll: list = lz
        ci: int = 1
        cj: int = 0
        cl: int = 2
    elif magnetic_component == 2:
        """Z Componente"""
        mli: list = ly
        mlj: list = lz
        mll: list = lx
        ci: int = 1
        cj: int = 2
        cl: int = 0
    elif magnetic_component not in [0, 1, 2]:
        raise ValueError(f"***Error\n\n Component not exist: {magnetic_component}")

    for i in range(total_nprim):

        for j in range(i, total_nprim):

            temp: float = 0.0
            if magnetic_component == 1:
                sxx: float = 0
                sxx = hermite_coefficient(
                                            lx[i],
                                            lx[j],
                                            0,
                                            coord[center[i]][0] - coord[center[j]][0],
                                            exp[i],
                                            exp[j],
                                        )
                syy = hermite_coefficient(
                            ly[i],
                            ly[j],
                            0,
                            coord[center[i]][1] - coord[center[j]][1],
                            exp[i],
                            exp[j],
                        )
                szz = hermite_coefficient(
                            lz[i],
                            lz[j],
                            0,
                            coord[center[i]][2] - coord[center[j]][2],
                            exp[i],
                            exp[j],
                        )
                d2phi: float = 0.0
                for sij, c, ml in zip([syy*szz, sxx*syy], [0,2], [lx,lz]):
                    d2phi += sij*(
                                4.0
                                * exp[j]
                                * exp[j]
                                * (
                                    hermite_coefficient(
                                        ml[i],
                                        ml[j] + 2,
                                        0,
                                        coord[center[i]][c] - coord[center[j]][c],
                                        exp[i],
                                        exp[j],
                                    )
                                )
                                - 2.0 * exp[j] * (2.0 * ml[j] + 1.0) 
                                * hermite_coefficient(
                                                        ml[i],
                                                        ml[j],
                                                        0,
                                                        coord[center[i]][c] - coord[center[j]][c],
                                                        exp[i],
                                                        exp[j],
                                                    ) 
                                + ml[j]
                                * (ml[j] - 1.0)
                                * (
                                    hermite_coefficient(
                                        ml[i],
                                        ml[j] - 2,
                                        0,
                                        coord[center[i]][c] - coord[center[j]][c],
                                        exp[i],
                                        exp[j],
                                    )
                                )
                            )
                temp = d2phi
            else:
                sll: float = 0.0
                sll = hermite_coefficient(
                                            mll[i],
                                            mll[j],
                                            0,
                                            coord[center[i]][cl] - coord[center[j]][cl],
                                            exp[i],
                                            exp[j],
                                        )
                didjphi: float = 1.0
                for c, ml in zip([ci,cj],[mli,mlj]):
                    didjphi *= (2.0 * exp[j] * hermite_coefficient(
                                ml[i],
                                ml[j] + 1,
                                0,
                                coord[center[i]][c]-coord[center[j]][c],
                                exp[i],
                                exp[j],
                            ) - ml[j] * hermite_coefficient(
                                ml[i],
                                ml[j] - 1,
                                0,
                                coord[center[i]][c]-coord[center[j]][c],
                                exp[i],
                                exp[j],
                            ))
                temp = -sll*didjphi

            curllgyp[count] = (
                normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * temp 
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )
            count += 1

    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Curl Angular Momentum Atomic Integrals on Y by Lineal Momentum on {magnetic_component} ",
            delta_time=(time() - start),
        )

    return curllgyp

if __name__ == "__main__":
    # STO-2G
    s = curllgyp(
        coord=[[0.0, 0.0, -0.545857052],[0.0, 0.0, 2.309057052]],
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
        magnetic_component=2,
        output=9,
        dalton_normalization=False,
        driver_time=None,
    )

    print("\n curllygp : ", s, "\n", len(s))
