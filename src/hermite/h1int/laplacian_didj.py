from lib1h import *

############# Calculate the potential one body integrals ########################
def laplacian_didj(
    coord, component, exp, center, lx, ly, lz, output, dalton_normalization, driver_time
):
    """_summary_

    Laplacian dd/didj

    Args:
        coord (list): list 2d with coordinates of the atoms
        component (int): Component to take of double derivative
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
        output (int): Output level for integral calculation
        dalton_normalization (bool): it is used the dalton normalization formule
        drive_time (drv_object): Object to manage the time

    Return:
        ddi (array): array 1d with atomic integrals
    """

    start = time()
    # Primitive total in the cluster
    total_nprim = len(exp)
    laplap = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]


    #% d4i or d2xd2y
    if component == 0 or component == 3 or component == 9 or component == 15: 
        # d4x or d2xd2y or d2xdydz or d3xdy
        mla: list = lx
        mlb: list = ly
        mlc: list = lz
        a: int = 0
        b: int = 1
        c: int = 2
    elif component == 1 or component == 5 or component == 11 or component == 17:
        # d4y or d2yd2x or d2ydxdz or d3ydx
        mla: list = ly
        mlb: list = lx
        mlc: list = lz
        a: int = 1
        b: int = 0
        c: int = 2
    elif component == 2 or component == 8 or component == 14 or component == 20:
        # d4z or d2zd2y or d2zdydx or d3zdy
        mla: list = lz
        mlb: list = ly
        mlc: list = lx
        a: int = 2
        b: int = 1
        c: int = 0
    elif component == 4 or component == 10 or component == 16:
        # d2xd2z or d2xdzdy or d3xdz
        mla: list = lx
        mlb: list = lz
        mlc: list = ly
        a: int = 0
        b: int = 2
        c: int = 1
    elif component == 6 or component == 12 or component == 18:
        # d2yd2z or d2ydzdx or d3ydz
        mla: list = ly
        mlb: list = lz
        mlc: list = lx
        a: int = 1
        b: int = 2
        c: int = 0
    elif component == 7 or component == 13 or component == 19:
        # d2zd2x or d3zdx
        mla: list = lz
        mlb: list = lx
        mlc: list = ly
        a: int = 2
        b: int = 0
        c: int = 1
    else:
        raise  ValueError(f"***Error\n\n This second derivative not exist: {component}")

    count = 0
    for i in range(total_nprim):

        for j in range(i, total_nprim):
            # Horizontal Reccurence
            # <phi|dxx|phi> = (dxxSij^0)Skl^0Smn^0 =
            # Dij^2Skl^0Smn^0 =
            # (4b^2Sij+2^0 - 2b(2j + 1)Sij^0 + j(j-1)Sij-2^0)Skl^0Smn^0 =
            # (4b^2E0^ij+2 - 2b(2j + 1)E0^ij + j(j-1)E0^ij-2 )E0^klE0^mn

            if component <= 9 or component >= 15:
                smn = hermite_coefficient(
                    mlc[i],
                    mlc[j],
                    0,
                    coord[center[i]][c] - coord[center[j]][c],
                    exp[i],
                    exp[j],
                )

            if component  >= 3 and component <= 14:  
                #% d2id2j or d2idjdk
                sij = (4.0 * exp[j] * exp[j] 
                        * hermite_coefficient(
                        mla[i],
                        mla[j] + 2,
                        0,
                        coord[center[i]][a] - coord[center[j]][a],
                        exp[i],
                        exp[j],
                        ) 
                        - 2.0 * exp[j] * (2.0 * mla[j] + 1.0) 
                        * hermite_coefficient(
                        mla[i],
                        mla[j],
                        0,
                        coord[center[i]][a] - coord[center[j]][a],
                        exp[i],
                        exp[j],
                        )
                        + mla[j] * (mla[j] - 1.0) 
                        * hermite_coefficient(
                        mla[i],
                        mla[j] - 2,
                        0,
                        coord[center[i]][a] - coord[center[j]][a],
                        exp[i],
                        exp[j],
                        )
                    )

            if component >= 9: 
                #% d2idjdk or d3idj
                skl = (2.0 * exp[j] 
                        * hermite_coefficient(
                                mlb[i],
                                mlb[j] + 1,
                                0,
                                coord[center[i]][b] - coord[center[j]][b],
                                exp[i],
                                exp[j],
                                ) 
                        - mlb[j]
                        * hermite_coefficient(
                            mlb[i],
                            mlb[j] - 1,
                            0,
                            coord[center[i]][b] - coord[center[j]][b],
                            exp[i],
                            exp[j],
                        )
                    )

            if component >= 0 and component <= 2:  
                #% dx^4, dy^4, dz^4
                skl = hermite_coefficient(
                    mlb[i],
                    mlb[j],
                    0,
                    coord[center[i]][b] - coord[center[j]][b],
                    exp[i],
                    exp[j],
                )
                sij = (
                        16.0 * exp[j] * exp[j] * exp[j] * exp[j]
                        * hermite_coefficient(
                                mla[i],
                                mla[j] + 4,
                                0,
                                coord[center[i]][component] - coord[center[j]][component],
                                exp[i],
                                exp[j],
                            )
                        - 16.0 * exp[j] * exp[j] * exp[j] * (2.0 * mla[j] + 3.0)
                        * hermite_coefficient(
                            mla[i],
                            mla[j] + 2,
                            0,
                            coord[center[i]][component] - coord[center[j]][component],
                            exp[i],
                            exp[j],
                        )
                        + 12.0 * exp[j] * exp[j] * (2.0 * mla[j] * mla[j] + 2.0 * mla[j] + 1.0)
                        * hermite_coefficient(
                                mla[i],
                                mla[j],
                                0,
                                coord[center[i]][component] - coord[center[j]][component],
                                exp[i],
                                exp[j],
                            )
                        - 4.0 * exp[j] * mla[j] * (mla[j] - 1.0) * (2.0 * mla[j] - 1.0)
                        * hermite_coefficient(
                                mla[i],
                                mla[j] - 2,
                                0,
                                coord[center[i]][component] - coord[center[j]][component],
                                exp[i],
                                exp[j],
                            )
                        + mla[j] * (mla[j] -1.0) * (mla[j] - 2.0) * (mla[j] - 3.0)
                        * hermite_coefficient(
                                mla[i],
                                mla[j] - 4,
                                0,
                                coord[center[i]][component] - coord[center[j]][component],
                                exp[i],
                                exp[j],
                            )
                        )
            elif component  >= 3 and component <= 8: 
                #% d2id2j
                skl = (4.0 * exp[j] * exp[j] 
                        * hermite_coefficient(
                        mlb[i],
                        mlb[j] + 2,
                        0,
                        coord[center[i]][b] - coord[center[j]][b],
                        exp[i],
                        exp[j],
                        ) 
                        - 2.0 * exp[j] * (2.0 * mlb[j] + 1.0) 
                        * hermite_coefficient(
                        mlb[i],
                        mlb[j],
                        0,
                        coord[center[i]][b] - coord[center[j]][b],
                        exp[i],
                        exp[j],
                        )
                        + mlb[j] * (mlb[j] - 1.0) 
                        * hermite_coefficient(
                        mlb[i],
                        mlb[j] - 2,
                        0,
                        coord[center[i]][b] - coord[center[j]][b],
                        exp[i],
                        exp[j],
                        )
                    )
            elif component >= 9 and component <= 14:
                #% d2idjdk
                smn = (2.0 * exp[j] 
                        * hermite_coefficient(
                                mlc[i],
                                mlc[j] + 1,
                                0,
                                coord[center[i]][c] - coord[center[j]][c],
                                exp[i],
                                exp[j],
                            ) 
                        - mlc[j]
                        * hermite_coefficient(
                            mlc[i],
                            mlc[j] - 1,
                            0,
                            coord[center[i]][c] - coord[center[j]][c],
                            exp[i],
                            exp[j],
                        )
                    )
            else:
                #% d3idj
                sij = (
                        8.0 * exp[j] * exp[j] * exp[j]
                        * hermite_coefficient(
                                mla[i],
                                mla[j] + 3,
                                0,
                                coord[center[i]][a] - coord[center[j]][a],
                                exp[i],
                                exp[j],
                            )
                        - 12.0 * exp[j] * exp[j] * (mla[j] + 1.0)
                        * hermite_coefficient(
                                mla[i],
                                mla[j] + 1,
                                0,
                                coord[center[i]][a] - coord[center[j]][a],
                                exp[i],
                                exp[j],
                            )
                        + 6.0 * exp[j] * mla[j] * mla[j]
                        * hermite_coefficient(
                                mla[i],
                                mla[j] - 1,
                                0,
                                coord[center[i]][a] - coord[center[j]][a],
                                exp[i],
                                exp[j],
                            )
                        - mla[i] * (mla[i] - 1.0) * (mla[i] - 2.0)
                        * hermite_coefficient(
                                mla[i],
                                mla[j] - 3,
                                0,
                                coord[center[i]][a] - coord[center[j]][a],
                                exp[i],
                                exp[j],
                            )
                      )

            laplap[count] = (
                  normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * (sij * skl * smn)
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )
            count += 1

    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Laplacian component {component}", delta_time=(time() - start)
        )

    return laplap
