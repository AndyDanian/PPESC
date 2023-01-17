from lib1h import *

############# Calculate the potential one body integrals ########################
def laplacian(
    coord, component, exp, center, lx, ly, lz, output, dalton_normalization, driver_time
):
    """_summary_

    Potential integrals

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
    ddi = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    if component == 0 or component == 3:
        mla: list = lx
        mlb: list = ly
        mlc: list = lz
        a: int = 0
        b: int = 1
        c: int = 2
    elif component == 1:
        mla: list = ly
        mlb: list = lx
        mlc: list = lz
        b: int = 0
        c: int = 2
    elif component == 2:
        mla: list = lz
        mlb: list = ly
        mlc: list = lx
        b: int = 1
        c: int = 0
    elif component == 4:
        mla: list = lx
        mlb: list = lz
        mlc: list = ly
        a: int = 0
        b: int = 2
        c: int = 1
    elif component == 5:
        mla: list = ly
        mlb: list = lz
        mlc: list = lx
        a: int = 1
        b: int = 2
        c: int = 0

    if component >= 0 and component <= 2:
        FSIGN: float = 1.0
    else:
        FSIGN: float = 1.0

    count = 0
    for i in range(total_nprim):

        for j in range(i, total_nprim):
            # Horizontal Reccurence
            # <phi|dxx|phi> = (dxxSij^0)Skl^0Smn^0 =
            # Dij^2Skl^0Smn^0 =
            # (4b^2Sij+2^0 - 2b(2j + 1)Sij^0 + j(j-1)Sij-2^0)Skl^0Smn^0 =
            # (4b^2E0^ij+2 - 2b(2j + 1)E0^ij + j(j-1)E0^ij-2 )E0^klE0^mn

            smn = hermite_coefficient(
                mlc[i],
                mlc[j],
                0,
                coord[center[i]][c] - coord[center[j]][c],
                exp[i],
                exp[j],
            )

            if component >= 0 and component <= 2:  # dx^2, dy^2, dz^2
                skl = hermite_coefficient(
                    mlb[i],
                    mlb[j],
                    0,
                    coord[center[i]][b] - coord[center[j]][b],
                    exp[i],
                    exp[j],
                )
                sij = (
                    4.0
                    * exp[j]
                    * exp[j]
                    * (
                        hermite_coefficient(
                            mla[i],
                            mla[j] + 2,
                            0,
                            coord[center[i]][component] - coord[center[j]][component],
                            exp[i],
                            exp[j],
                        )
                    )
                    - 2.0
                    * exp[j]
                    * (2.0 * mla[j] + 1.0)
                    * hermite_coefficient(
                        mla[i],
                        mla[j],
                        0,
                        coord[center[i]][component] - coord[center[j]][component],
                        exp[i],
                        exp[j],
                    )
                    + mla[j]
                    * (mla[j] - 1.0)
                    * (
                        hermite_coefficient(
                            mla[i],
                            mla[j] - 2,
                            0,
                            coord[center[i]][component] - coord[center[j]][component],
                            exp[i],
                            exp[j],
                        )
                    )
                )
            else:  # dxy, dxz, dyz
                # Horizontal Reccurence
                # <phi|dxdy|phi> = (dxSij^0)(dySkl^0)Smn^0 =
                # Dij^1Dkl^1Smn^0=
                # [dx(x-Ax)^ie^-a(x-Ax)^2][(x-Bx)^je^-b(x-Bx)^2]
                # [(y-Ay)^ke^-a(y-Ay)^2][dy(y-By)^le^-b(y-By)^2]Smn^0 =
                # [{i(x-Ax)^i-1 - 2a(x-Ax)^i+1}e^-a(x-Ax)^2][(x-Bx)^je^-b(x-Bx)^2]
                # [(y-Ay)^ke^-a(y-Ay)^2][{l(y-By)^l-1 - 2b(y-By)^l+1}e^-b(y-By)^2]Smn^0 =
                # (iSi-1j^0 - 2aSi+1j^0)(lSkl-1^0 - 2bSkl+1^0)Smn^0 =
                # (iE0^i-1j - 2aE0^i+1)(lE0^kl-1 - 2bE0^kl+1)E0^mn
                sij = 2.0 * exp[i] * hermite_coefficient(
                    mla[i] + 1,
                    mla[j],
                    0,
                    coord[center[i]][a] - coord[center[j]][a],
                    exp[i],
                    exp[j],
                ) - (
                    mla[i]
                    * (
                        hermite_coefficient(
                            mla[i] - 1,
                            mla[j],
                            0,
                            coord[center[i]][a] - coord[center[j]][a],
                            exp[i],
                            exp[j],
                        )
                    )
                )
                skl = 2.0 * exp[j] * hermite_coefficient(
                    mlb[i],
                    mlb[j] + 1,
                    0,
                    coord[center[i]][b] - coord[center[j]][b],
                    exp[i],
                    exp[j],
                ) - (
                    mlb[j]
                    * (
                        hermite_coefficient(
                            mlb[i],
                            mlb[j] - 1,
                            0,
                            coord[center[i]][b] - coord[center[j]][b],
                            exp[i],
                            exp[j],
                        )
                    )
                )

            ddi[count] = (
                FSIGN
                * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * (sij * skl * smn)
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )
            count += 1

    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Laplacian component {component}", delta_time=(time() - start)
        )

    return ddi
