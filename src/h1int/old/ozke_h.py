from libh import *

start = time()

# 6-311++G**
exp = [
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
]

n = [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1]
lx = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]
ly = [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
lz = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1]

center = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1]
total_nprim = 18

coord = [[0.0, 0.0, 0.0586476414], [0.0, 0.0, 1.4045523587]]

intLx = np.zeros((total_nprim, total_nprim), dtype=float)
intLy = np.zeros((total_nprim, total_nprim), dtype=float)
intLz = np.zeros((total_nprim, total_nprim), dtype=float)
output = 11

# !Note: Rgaugeo is put on in the origen
Rg = [0.0, 0.0, 1.4045523587]
for i in range(total_nprim):

    for j in range(i, total_nprim):

        Px = (
            exp[i] * coord[center[i]][0]
            + exp[j] * coord[center[j]][0]
        )
        Px = Px / (exp[i] + exp[j])

        Py = (
            exp[i] * coord[center[i]][1]
            + exp[j] * coord[center[j]][1]
        )
        Py = Py / (exp[i] + exp[j])

        Pz = (
            exp[i] * coord[center[i]][2]
            + exp[j] * coord[center[j]][2]
        )
        Pz = Pz / (exp[i] + exp[j])

        Xpg = Px - Rg[0]
        Ypg = Py - Rg[1]
        Zpg = Pz - Rg[2]

        # E_0^ij
        e0ij = E(
            lx[i],
            lx[j],
            0,
            coord[center[i]][0] - coord[center[j]][0],
            exp[i],
            exp[j],
        )

        # E_0^kl
        e0kl = E(
            ly[i],
            ly[j],
            0,
            coord[center[i]][1] - coord[center[j]][1],
            exp[i],
            exp[j],
        )

        # E_0^mn
        e0mn = E(
            lz[i],
            lz[j],
            0,
            coord[center[i]][2] - coord[center[j]][2],
            exp[i],
            exp[j],
        )

        # E_1^ij
        e1ij = E(
            lx[i],
            lx[j],
            1,
            coord[center[i]][0] - coord[center[j]][0],
            exp[i],
            exp[j],
        )

        # E_1^kl
        e1kl = E(
            ly[i],
            ly[j],
            1,
            coord[center[i]][1] - coord[center[j]][1],
            exp[i],
            exp[j],
        )

        # E_1^mn
        e1mn = E(
            lz[i],
            lz[j],
            1,
            coord[center[i]][2] - coord[center[j]][2],
            exp[i],
            exp[j],
        )

        # D_ij^1
        px = 2.0 * exp[j] * E(
            lx[i],
            lx[j] + 1,
            0,
            coord[center[i]][0] - coord[center[j]][0],
            exp[i],
            exp[j],
        ) - lx[j] * E(
            lx[i],
            lx[j] - 1,
            0,
            coord[center[i]][0] - coord[center[j]][0],
            exp[i],
            exp[j],
        )

        # D_kl^1
        # ! izquierda
        yp = 2.0 * exp[i] * E(
            ly[i] + 1,
            ly[j],
            0,
            coord[center[i]][1] - coord[center[j]][1],
            exp[i],
            exp[j],
        ) - ly[i] * E(
            ly[i] - 1,
            ly[j],
            0,
            coord[center[i]][1] - coord[center[j]][1],
            exp[i],
            exp[j],
        )
        # ! derecha
        py = 2.0 * exp[j] * E(
            ly[i],
            ly[j] + 1,
            0,
            coord[center[i]][1] - coord[center[j]][1],
            exp[i],
            exp[j],
        ) - ly[j] * E(
            ly[i],
            ly[j] - 1,
            0,
            coord[center[i]][1] - coord[center[j]][1],
            exp[i],
            exp[j],
        )

        # D_mn^1
        # ! izquierda
        zp = 2.0 * exp[i] * E(
            lz[i] + 1,
            lz[j],
            0,
            coord[center[i]][2] - coord[center[j]][2],
            exp[i],
            exp[j],
        ) - lz[i] * E(
            lz[i] - 1,
            lz[j],
            0,
            coord[center[i]][2] - coord[center[j]][2],
            exp[i],
            exp[j],
        )
        # ! derecha
        pz = 2.0 * exp[j] * E(
            lz[i],
            lz[j] + 1,
            0,
            coord[center[i]][2] - coord[center[j]][2],
            exp[i],
            exp[j],
        ) - lz[j] * E(
            lz[i],
            lz[j] - 1,
            0,
            coord[center[i]][2] - coord[center[j]][2],
            exp[i],
            exp[j],
        )

        # ! Terms of dxxLx
        xxd = (
            4.0
            * exp[i]
            * exp[i]
            * (
                E(
                    lx[i] + 2,
                    lx[j],
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
            )
            - 2.0 * exp[i] * (2.0 * lx[i] + 1.0) * e0ij
            + lx[i]
            * (lx[i] - 1.0)
            * (
                E(
                    lx[i] - 2,
                    lx[j],
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
            )
        )

        dxxlx = xxd * ((e1kl + Ypg * e0kl) * pz - (e1mn + Zpg * e0mn) * py)

        # ! Terms of dyyLx
        ########### Dipole
        # E_1^kl + Ypc*E_0^kl
        d_skl_1 = 2.0 * exp[i] * (2.0 * ly[i] + 1) * (e1kl + Ypg * e0kl)
        # E_1^k-2l + Ypc*E_0^k-2l
        d_skm2l_1 = (
            ly[i]
            * (ly[i] - 1)
            * (
                E(
                    ly[i] - 2,
                    ly[j],
                    1,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
                + Ypg
                * E(
                    ly[i] - 2,
                    ly[j],
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
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
                E(
                    ly[i] + 2,
                    ly[j],
                    1,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
                + Ypg
                * E(
                    ly[i] + 2,
                    ly[j],
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
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
            * (2.0 * ly[i] + 1.0)
            * (
                2.0
                * exp[j]
                * E(
                    ly[i],
                    ly[j] + 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
                - ly[j]
                * E(
                    ly[i],
                    ly[j] - 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
            )
        )
        # E_0^k-2l+1 - l*E_0^k-2l-1
        d_dkm2l_1 = (
            ly[i]
            * (ly[i] - 1.0)
            * (
                2.0
                * exp[j]
                * E(
                    ly[i] - 2,
                    ly[j] + 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
                - ly[j]
                * E(
                    ly[i] - 2,
                    ly[j] - 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
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
                * E(
                    ly[i] + 2,
                    ly[j] + 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
                - ly[j]
                * E(
                    ly[i] + 2,
                    ly[j] - 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
            )
        )

        dyylx = (
            pz * (d_skm2l_1 - d_skl_1 + d_skt2l_1)
            + (e1mn + Zpg * e0mn) * (d_dkl_1 - d_dkm2l_1 - d_dkt2l_1)
        ) * e0ij

        # ! Terms of dzzLx
        # E_1^mn + Zpc*E_0^mn
        d_smn_1 = 2.0 * exp[i] * (2.0 * lz[i] + 1) * (e1mn + Zpg * e0mn)
        # E_1^m-2n + Zpc*E_0^m-2n
        d_smm2n_1 = (
            lz[i]
            * (lz[i] - 1)
            * (
                E(
                    lz[i] - 2,
                    lz[j],
                    1,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
                + Zpg
                * E(
                    lz[i] - 2,
                    lz[j],
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
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
                E(
                    lz[i] + 2,
                    lz[j],
                    1,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
                + Zpg
                * E(
                    lz[i] + 2,
                    lz[j],
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
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
            * (2.0 * lz[i] + 1.0)
            * (
                2.0
                * exp[j]
                * E(
                    lz[i],
                    lz[j] + 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
                - lz[j]
                * E(
                    lz[i],
                    lz[j] - 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
            )
        )
        # E_0^k-2l+1 - l*E_0^k-2l-1
        d_dmm2n_1 = (
            lz[i]
            * (lz[i] - 1.0)
            * (
                2.0
                * exp[j]
                * E(
                    lz[i] - 2,
                    lz[j] + 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
                - lz[j]
                * E(
                    lz[i] - 2,
                    lz[j] - 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
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
                * E(
                    lz[i] + 2,
                    lz[j] + 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
                - lz[j]
                * E(
                    lz[i] + 2,
                    lz[j] - 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
            )
        )

        dzzlx = (
            py * (d_smn_1 - d_smm2n_1 - d_smt2n_1)
            + (e1kl + Ypg * e0kl) * (d_dmm2n_1 - d_dmn_1 + d_dmt2n_1)
        ) * e0ij

        # * nabla Real{Lx} + Real{Lx} nabla

        intLx[i, j] = (
            Norm[lx[i] + ly[i] + lz[i]](exp[i])
            * Norm[lx[j] + ly[j] + lz[j]](exp[j])
            * 0.25
            * (dxxlx + dyylx + dzzlx)
            * np.power(np.pi / (exp[i] + exp[j]), 1.5)
        )
        intLx[j, i] = -1.0 * intLx[i, j]

        ###########################
        # ! Terms of dyyLy
        yyd = (
            4.0
            * exp[i]
            * exp[i]
            * (
                E(
                    ly[i] + 2,
                    ly[j],
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
            )
            - 2.0 * exp[i] * (2.0 * ly[i] + 1.0) * e0kl
            + ly[i]
            * (ly[i] - 1.0)
            * (
                E(
                    ly[i] - 2,
                    ly[j],
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
            )
        )

        dyyly = yyd * ((e1mn + Zpg * e0mn) * px - (e1ij + Xpg * e0ij) * pz)

        # ! Terms of dzzLy
        ########### Dipole
        # E_1^kl + Ypc*E_0^kl
        d_skl_1 = 2.0 * exp[i] * (2.0 * lz[i] + 1) * (e1mn + Zpg * e0mn)
        # E_1^k-2l + Ypc*E_0^k-2l
        d_skm2l_1 = (
            lz[i]
            * (lz[i] - 1)
            * (
                E(
                    lz[i] - 2,
                    lz[j],
                    1,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
                + Zpg
                * E(
                    lz[i] - 2,
                    lz[j],
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
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
                E(
                    lz[i] + 2,
                    lz[j],
                    1,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
                + Zpg
                * E(
                    lz[i] + 2,
                    lz[j],
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
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
            * (2.0 * lz[i] + 1.0)
            * (
                2.0
                * exp[j]
                * E(
                    lz[i],
                    lz[j] + 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
                - lz[j]
                * E(
                    lz[i],
                    lz[j] - 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
            )
        )
        # E_0^k-2l+1 - l*E_0^k-2l-1
        d_dkm2l_1 = (
            lz[i]
            * (lz[i] - 1.0)
            * (
                2.0
                * exp[j]
                * E(
                    lz[i] - 2,
                    lz[j] + 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
                - lz[j]
                * E(
                    lz[i] - 2,
                    lz[j] - 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
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
                * E(
                    lz[i] + 2,
                    lz[j] + 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
                - lz[j]
                * E(
                    lz[i] + 2,
                    lz[j] - 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
            )
        )

        dzzly = (
            px * (d_skm2l_1 - d_skl_1 + d_skt2l_1)
            + (e1ij + Xpg * e0ij) * (d_dkl_1 - d_dkm2l_1 - d_dkt2l_1)
        ) * e0kl
        # ! Terms of dxxLy
        # E_1^mn + Zpc*E_0^mn
        d_smn_1 = 2.0 * exp[i] * (2.0 * lx[i] + 1) * (e1ij + Xpg * e0ij)
        # E_1^m-2n + Zpc*E_0^m-2n
        d_smm2n_1 = (
            lx[i]
            * (lx[i] - 1)
            * (
                E(
                    lx[i] - 2,
                    lx[j],
                    1,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
                + Xpg
                * E(
                    lx[i] - 2,
                    lx[j],
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
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
                E(
                    lx[i] + 2,
                    lx[j],
                    1,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
                + Xpg
                * E(
                    lx[i] + 2,
                    lx[j],
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
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
            * (2.0 * lx[i] + 1.0)
            * (
                2.0
                * exp[j]
                * E(
                    lx[i],
                    lx[j] + 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
                - lx[j]
                * E(
                    lx[i],
                    lx[j] - 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
            )
        )
        # E_0^k-2l+1 - l*E_0^k-2l-1
        d_dmm2n_1 = (
            lx[i]
            * (lx[i] - 1.0)
            * (
                2.0
                * exp[j]
                * E(
                    lx[i] - 2,
                    lx[j] + 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
                - lx[j]
                * E(
                    lx[i] - 2,
                    lx[j] - 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
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
                * E(
                    lx[i] + 2,
                    lx[j] + 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
                - lx[j]
                * E(
                    lx[i] + 2,
                    lx[j] - 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
            )
        )

        dxxly = (
            pz * (d_smn_1 - d_smm2n_1 - d_smt2n_1)
            + (e1mn + Zpg * e0mn) * (d_dmm2n_1 - d_dmn_1 + d_dmt2n_1)
        ) * e0kl
        # * nabla Real{Ly} + Real{Ly} nabla

        intLy[i, j] = (
            Norm[lx[i] + ly[i] + lz[i]](exp[i])
            * Norm[lx[j] + ly[j] + lz[j]](exp[j])
            * 0.25
            * (dxxly + dyyly + dzzly)
            * np.power(np.pi / (exp[i] + exp[j]), 1.5)
        )
        intLy[j, i] = -1.0 * intLy[i, j]

        ###########################
        # ! Terms of dzzLz
        zzd = (
            4.0
            * exp[i]
            * exp[i]
            * (
                E(
                    lz[i] + 2,
                    lz[j],
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
            )
            - 2.0 * exp[i] * (2.0 * lz[i] + 1.0) * e0mn
            + lz[i]
            * (lz[i] - 1.0)
            * (
                E(
                    lz[i] - 2,
                    lz[j],
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp[i],
                    exp[j],
                )
            )
        )

        dzzlz = zzd * ((e1ij + Xpg * e0ij) * py - (e1kl + Ypg * e0kl) * px)

        # ! Terms of dxxLz
        ########### Dipole

        d_skl_1 = 2.0 * exp[i] * (2.0 * lx[i] + 1) * (e1ij + Xpg * e0ij)

        d_skm2l_1 = (
            lx[i]
            * (lx[i] - 1)
            * (
                E(
                    lx[i] - 2,
                    lx[j],
                    1,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
                + Xpg
                * E(
                    lx[i] - 2,
                    lx[j],
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
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
                E(
                    lx[i] + 2,
                    lx[j],
                    1,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
                + Xpg
                * E(
                    lx[i] + 2,
                    lx[j],
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
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
            * (2.0 * lx[i] + 1.0)
            * (
                2.0
                * exp[j]
                * E(
                    lx[i],
                    lx[j] + 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
                - lx[j]
                * E(
                    lx[i],
                    lx[j] - 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
            )
        )

        d_dkm2l_1 = (
            lx[i]
            * (lx[i] - 1.0)
            * (
                2.0
                * exp[j]
                * E(
                    lx[i] - 2,
                    lx[j] + 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
                - lx[j]
                * E(
                    lx[i] - 2,
                    lx[j] - 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
            )
        )

        d_dkt2l_1 = (
            4.0
            * exp[i]
            * exp[i]
            * (
                2.0
                * exp[j]
                * E(
                    lx[i] + 2,
                    lx[j] + 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
                - lx[j]
                * E(
                    lx[i] + 2,
                    lx[j] - 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp[i],
                    exp[j],
                )
            )
        )

        dxxlz = (
            py * (d_skm2l_1 - d_skl_1 + d_skt2l_1)
            + (e1kl + Ypg * e0kl) * (d_dkl_1 - d_dkm2l_1 - d_dkt2l_1)
        ) * e0mn

        # ! Terms of dyyLz

        d_smn_1 = 2.0 * exp[i] * (2.0 * ly[i] + 1) * (e1kl + Ypg * e0kl)

        d_smm2n_1 = (
            ly[i]
            * (ly[i] - 1)
            * (
                E(
                    ly[i] - 2,
                    ly[j],
                    1,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
                + Ypg
                * E(
                    ly[i] - 2,
                    ly[j],
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
            )
        )

        d_smt2n_1 = (
            4.0
            * exp[i]
            * exp[i]
            * (
                E(
                    ly[i] + 2,
                    ly[j],
                    1,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
                + Ypg
                * E(
                    ly[i] + 2,
                    ly[j],
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
            )
        )
        ####################### Derivatives

        d_dmn_1 = (
            2.0
            * exp[i]
            * (2.0 * ly[i] + 1.0)
            * (
                2.0
                * exp[j]
                * E(
                    ly[i],
                    ly[j] + 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
                - ly[j]
                * E(
                    ly[i],
                    ly[j] - 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
            )
        )

        d_dmm2n_1 = (
            ly[i]
            * (ly[i] - 1.0)
            * (
                2.0
                * exp[j]
                * E(
                    ly[i] - 2,
                    ly[j] + 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
                - ly[j]
                * E(
                    ly[i] - 2,
                    ly[j] - 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
            )
        )

        d_dmt2n_1 = (
            4.0
            * exp[i]
            * exp[i]
            * (
                2.0
                * exp[j]
                * E(
                    ly[i] + 2,
                    ly[j] + 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
                - ly[j]
                * E(
                    ly[i] + 2,
                    ly[j] - 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp[i],
                    exp[j],
                )
            )
        )

        dyylz = (
            px * (d_smn_1 - d_smm2n_1 - d_smt2n_1)
            + (e1ij + Xpg * e0ij) * (d_dmm2n_1 - d_dmn_1 + d_dmt2n_1)
        ) * e0mn
        # * nabla Real{Lz} + Real{Lz} nabla

        intLz[i, j] = (
            Norm[lx[i] + ly[i] + lz[i]](exp[i])
            * Norm[lx[j] + ly[j] + lz[j]](exp[j])
            * 0.25
            * (dxxlz + dyylz + dzzlz)
            * np.power(np.pi / (exp[i] + exp[j]), 1.5)
        )
        intLz[j, i] = -1.0 * intLz[i, j]

        if (
            output > 10
            and np.abs(intLy[i, j]) > 1e-3  # and j == 7 and i == 0
        ):  # :
            print(
                "int [",
                i + 1,
                ",",
                j + 1,
                "] : ",
                dxxly, dyyly, dzzly
            )

print(" time [s]: ", -start + time())
