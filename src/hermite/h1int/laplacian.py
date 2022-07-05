from lib1h import *

############# Calculate the potential one body integrals ########################
def laplacian(coord, component, exp, center, lx, ly, lz, output, dalton_normalization, driver_time):
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

    if component == 0:
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
        a: int = 1
        b: int = 0
        c: int = 2
    elif component == 2:
        mla: list = lz
        mlb: list = ly
        mlc: list = lx
        a: int = 2
        b: int = 1
        c: int = 0

    count = 0
    for i in range(total_nprim):

        for j in range(i, total_nprim):
            # Horizontal Reccurence
            # <phi|dxx|phi> = (dxxSij^0)Skl^0Smn^0 =
            # Dij^2Skl^0Smn^0 =
            # (4b^2Sij+2^0 - 2b(2j + 1)Sij^0 + j(j-1)Sij-2^0)Skl^0Smn^0 =
            # (4b^2E0^ij+2 - 2b(2j + 1)E0^ij + j(j-1)E0^ij-2 )E0^klE0^mn

            sij = E(
                mla[i],
                mla[j],
                0,
                coord[center[i]][a] - coord[center[j]][a],
                exp[i],
                exp[j],
            )

            skl = E(
                mlb[i],
                mlb[j],
                0,
                coord[center[i]][b] - coord[center[j]][b],
                exp[i],
                exp[j],
            )

            smn = E(
                mlc[i],
                mlc[j],
                0,
                coord[center[i]][c] - coord[center[j]][c],
                exp[i],
                exp[j],
            )

            ddis = (
                4.0
                * exp[j]
                * exp[j]
                * (
                    E(
                        mla[i],
                        mla[j] + 2,
                        0,
                        coord[center[i]][component] - coord[center[j]][component],
                        exp[i],
                        exp[j],
                    )
                )
                - 2.0 * exp[j] * (2.0 * mla[j] + 1.0) * sij
                + mla[j]
                * (mla[j] - 1.0)
                * (
                    E(
                        mla[i],
                        mla[j] - 2,
                        0,
                        coord[center[i]][component] - coord[center[j]][component],
                        exp[i],
                        exp[j],
                    )
                )
            )

            ddi[count] = (
                normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * (ddis * skl * smn)
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )
            count += 1

    if output > 10:
        driver_time.add_name_delta_time(name = f"Laplacian component {component}", delta_time = (time() - start))

    return ddi