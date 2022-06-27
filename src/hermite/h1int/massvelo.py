from lib1h import *

def massvelo(coord, exp, center, lx, ly, lz, output, dalton_normalization):
    """
    Spin dipolar atomic integrals, which is a tensor

    Args:
        coord (list): list 2d with coordinates of the atoms
        atom (int): atomic index
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
        output (int): Output level for integral calculation
        dalton_normalization (bool): it is used the dalton normalization formule

    Return:
        angmom (array): array 2d with atomic integrals
    """
    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    massvelo: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    count: int = 0

    SPEED_LIGHT: float = 137.0359998
    CONST_ALPHA2: float = 1.0 / (SPEED_LIGHT * SPEED_LIGHT)

    #! Note: It's neccesary to multiplicate by 5.0/3.0 to get DALTON's values
    #!       there is difference in the 6-th decimal
    CONST_MASSV: float = CONST_ALPHA2 / 8.0 * 5.0 / 3.0
    for i in range(total_nprim):

        for j in range(i, total_nprim):
            # Horizontal Reccurence
            # -(dxx<phi|)(dxx|phi>) =
            # (4a^2x^i+2 - 2a(2i + 1)x^i + i(i-1)x^i-2)y^kz^m
            # (4b^2x^j+2 - 2b(2j + 1)x^j + j(j-1)x^j-2)y^lz^n =
            # 16a^2b^2Si+2j+2^0 + i(i+1)j(j+1)Si-2j-2^0 +4ab(2i+1)(2j+1)Sij^0
            # + 4a^2j(j-1)Si+2j-2^0 + 4b^2i(i-1)Si-2j+2^0
            # - 8a^2b(2j+1)Si+2j^0 - 8b^2a(2i+1)Sij+2^0
            # - 2a(2i+1)j(j-1)Sij-2^0 - 2b(2j+1)i(i-1)Si-2j^0

            sij = E(
                lx[i],
                lx[j],
                0,
                coord[center[i]][0] - coord[center[j]][0],
                exp[i],
                exp[j],
            )

            skl = E(
                ly[i],
                ly[j],
                0,
                coord[center[i]][1] - coord[center[j]][1],
                exp[i],
                exp[j],
            )

            smn = E(
                lz[i],
                lz[j],
                0,
                coord[center[i]][2] - coord[center[j]][2],
                exp[i],
                exp[j],
            )

            # dxxdxx
            dxxsidxxsj = (
                16.0
                * exp[i]
                * exp[i]
                * exp[j]
                * exp[j]
                * (
                    E(
                        lx[i] + 2,
                        lx[j] + 2,
                        0,
                        coord[center[i]][0] - coord[center[j]][0],
                        exp[i],
                        exp[j],
                    )
                )
                + 4.0
                * exp[i]
                * exp[j]
                * (2.0 * lx[i] + 1.0)
                * (2.0 * lx[j] + 1.0)
                * sij
                + lx[i]
                * (lx[i] - 1.0)
                * lx[j]
                * (lx[j] - 1.0)
                * (
                    E(
                        lx[i] - 2,
                        lx[j] - 2,
                        0,
                        coord[center[i]][0] - coord[center[j]][0],
                        exp[i],
                        exp[j],
                    )
                )
                + 4.0
                * exp[i]
                * exp[i]
                * lx[j]
                * (lx[j] - 1.0)
                * E(
                    lx[i] + 2,
                    lx[j] - 2,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
                + 4.0
                * exp[j]
                * exp[j]
                * lx[i]
                * (lx[i] - 1.0)
                * E(
                    lx[i] - 2,
                    lx[j] + 2,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
                - 8.0
                * exp[i]
                * exp[i]
                * exp[j]
                * (2.0 * lx[j] + 1.0)
                * E(
                    lx[i] + 2,
                    lx[j],
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
                - 8.0
                * exp[j]
                * exp[j]
                * exp[i]
                * (2.0 * lx[i] + 1.0)
                * E(
                    lx[i],
                    lx[j] + 2,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
                - 2.0
                * exp[i]
                * (2.0 * lx[i] + 1)
                * lx[j]
                * (lx[j] - 1)
                * E(
                    lx[i],
                    lx[j] - 2,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
                - 2.0
                * exp[j]
                * (2.0 * lx[j] + 1)
                * lx[i]
                * (lx[i] - 1)
                * E(
                    lx[i] - 2,
                    lx[j],
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
            )

            # dyydyy,
            dyyskdyysl = (
                16.0
                * exp[i]
                * exp[i]
                * exp[j]
                * exp[j]
                * (
                    E(
                        ly[i] + 2,
                        ly[j] + 2,
                        0,
                        coord[center[i]][1] - coord[center[j]][1],
                        exp[i],
                        exp[j],
                    )
                )
                + 4.0
                * exp[i]
                * exp[j]
                * (2.0 * ly[i] + 1.0)
                * (2.0 * ly[j] + 1.0)
                * skl
                + ly[i]
                * (ly[i] - 1.0)
                * ly[j]
                * (ly[j] - 1.0)
                * (
                    E(
                        ly[i] - 2,
                        ly[j] - 2,
                        0,
                        coord[center[i]][1] - coord[center[j]][1],
                        exp[i],
                        exp[j],
                    )
                )
                + 4.0
                * exp[i]
                * exp[i]
                * ly[j]
                * (ly[j] - 1.0)
                * E(
                    ly[i] + 2,
                    ly[j] - 2,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
                + 4.0
                * exp[j]
                * exp[j]
                * ly[i]
                * (ly[i] - 1.0)
                * E(
                    ly[i] - 2,
                    ly[j] + 2,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
                - 8.0
                * exp[i]
                * exp[i]
                * exp[j]
                * (2.0 * ly[j] + 1.0)
                * E(
                    ly[i] + 2,
                    ly[j],
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
                - 8.0
                * exp[j]
                * exp[j]
                * exp[i]
                * (2.0 * ly[i] + 1.0)
                * E(
                    ly[i],
                    ly[j] + 2,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
                - 2.0
                * exp[i]
                * (2.0 * ly[i] + 1)
                * ly[j]
                * (ly[j] - 1)
                * E(
                    ly[i],
                    ly[j] - 2,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
                - 2.0
                * exp[j]
                * (2.0 * ly[j] + 1)
                * ly[i]
                * (ly[i] - 1)
                * E(
                    ly[i] - 2,
                    ly[j],
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
            )

            # dzzdzz,
            dzzsmdzzsn = (
                16.0
                * exp[i]
                * exp[i]
                * exp[j]
                * exp[j]
                * (
                    E(
                        lz[i] + 2,
                        lz[j] + 2,
                        0,
                        coord[center[i]][2] - coord[center[j]][2],
                        exp[i],
                        exp[j],
                    )
                )
                + 4.0
                * exp[i]
                * exp[j]
                * (2.0 * lz[i] + 1.0)
                * (2.0 * lz[j] + 1.0)
                * smn
                + lz[i]
                * (lz[i] - 1.0)
                * lz[j]
                * (lz[j] - 1.0)
                * (
                    E(
                        lz[i] - 2,
                        lz[j] - 2,
                        0,
                        coord[center[i]][2] - coord[center[j]][2],
                        exp[i],
                        exp[j],
                    )
                )
                + 4.0
                * exp[i]
                * exp[i]
                * lz[j]
                * (lz[j] - 1.0)
                * E(
                    lz[i] + 2,
                    lz[j] - 2,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
                + 4.0
                * exp[j]
                * exp[j]
                * lz[i]
                * (lz[i] - 1.0)
                * E(
                    lz[i] - 2,
                    lz[j] + 2,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
                - 8.0
                * exp[i]
                * exp[i]
                * exp[j]
                * (2.0 * lz[j] + 1.0)
                * E(
                    lz[i] + 2,
                    lz[j],
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
                - 8.0
                * exp[j]
                * exp[j]
                * exp[i]
                * (2.0 * lz[i] + 1.0)
                * E(
                    lz[i],
                    lz[j] + 2,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
                - 2.0
                * exp[i]
                * (2.0 * lz[i] + 1)
                * lz[j]
                * (lz[j] - 1)
                * E(
                    lz[i],
                    lz[j] - 2,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
                - 2.0
                * exp[j]
                * (2.0 * lz[j] + 1)
                * lz[i]
                * (lz[i] - 1)
                * E(
                    lz[i] - 2,
                    lz[j],
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
            )

            massvelo[count] = (
                -normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * CONST_MASSV
                * (
                    dxxsidxxsj * skl * smn
                    + sij * dyyskdyysl * smn
                    + sij * skl * dzzsmdzzsn
                )
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )
            count += 1

    if output > 0:
        print_time(
            name = f"Massvelo Atomic Integrals", delta_time = (time() - start)
        )

    return massvelo