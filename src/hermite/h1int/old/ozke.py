from numpy import exp
import numpy as np
import time
import phi

start = time.time()

# 6-311++G**
exp_array = [
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

Norm = {0: phi.NS, 1: phi.NP}
intLx = np.zeros((total_nprim, total_nprim), dtype=float)
intLy = np.zeros((total_nprim, total_nprim), dtype=float)
intLz = np.zeros((total_nprim, total_nprim), dtype=float)
output = 11

# !Note: Rgaugeo is put on in the origen
Rg = [0.0, 0.0, 1.4045523587]
for i in range(total_nprim):

    for j in range(i, total_nprim):

        Px = (
            exp_array[i] * coord[center[i]][0]
            + exp_array[j] * coord[center[j]][0]
        )
        Px = Px / (exp_array[i] + exp_array[j])

        Py = (
            exp_array[i] * coord[center[i]][1]
            + exp_array[j] * coord[center[j]][1]
        )
        Py = Py / (exp_array[i] + exp_array[j])

        Pz = (
            exp_array[i] * coord[center[i]][2]
            + exp_array[j] * coord[center[j]][2]
        )
        Pz = Pz / (exp_array[i] + exp_array[j])

        Xpg = Px - Rg[0]
        Ypg = Py - Rg[1]
        Zpg = Pz - Rg[2]

        # E_0^ij
        e0ij = phi.E(
            lx[i],
            lx[j],
            0,
            coord[center[i]][0] - coord[center[j]][0],
            exp_array[i],
            exp_array[j],
        )

        # E_0^kl
        e0kl = phi.E(
            ly[i],
            ly[j],
            0,
            coord[center[i]][1] - coord[center[j]][1],
            exp_array[i],
            exp_array[j],
        )

        # E_0^mn
        e0mn = phi.E(
            lz[i],
            lz[j],
            0,
            coord[center[i]][2] - coord[center[j]][2],
            exp_array[i],
            exp_array[j],
        )

        # E_1^ij
        e1ij = phi.E(
            lx[i],
            lx[j],
            1,
            coord[center[i]][0] - coord[center[j]][0],
            exp_array[i],
            exp_array[j],
        )

        # E_1^kl
        e1kl = phi.E(
            ly[i],
            ly[j],
            1,
            coord[center[i]][1] - coord[center[j]][1],
            exp_array[i],
            exp_array[j],
        )

        # E_1^mn
        e1mn = phi.E(
            lz[i],
            lz[j],
            1,
            coord[center[i]][2] - coord[center[j]][2],
            exp_array[i],
            exp_array[j],
        )

        # D_ij^1
        px = 2.0 * exp_array[j] * phi.E(
            lx[i],
            lx[j] + 1,
            0,
            coord[center[i]][0] - coord[center[j]][0],
            exp_array[i],
            exp_array[j],
        ) - lx[j] * phi.E(
            lx[i],
            lx[j] - 1,
            0,
            coord[center[i]][0] - coord[center[j]][0],
            exp_array[i],
            exp_array[j],
        )

        # D_kl^1
        # ! izquierda
        yp = 2.0 * exp_array[i] * phi.E(
            ly[i] + 1,
            ly[j],
            0,
            coord[center[i]][1] - coord[center[j]][1],
            exp_array[i],
            exp_array[j],
        ) - ly[i] * phi.E(
            ly[i] - 1,
            ly[j],
            0,
            coord[center[i]][1] - coord[center[j]][1],
            exp_array[i],
            exp_array[j],
        )
        # ! derecha
        py = 2.0 * exp_array[j] * phi.E(
            ly[i],
            ly[j] + 1,
            0,
            coord[center[i]][1] - coord[center[j]][1],
            exp_array[i],
            exp_array[j],
        ) - ly[j] * phi.E(
            ly[i],
            ly[j] - 1,
            0,
            coord[center[i]][1] - coord[center[j]][1],
            exp_array[i],
            exp_array[j],
        )

        # D_mn^1
        # ! izquierda
        zp = 2.0 * exp_array[i] * phi.E(
            lz[i] + 1,
            lz[j],
            0,
            coord[center[i]][2] - coord[center[j]][2],
            exp_array[i],
            exp_array[j],
        ) - lz[i] * phi.E(
            lz[i] - 1,
            lz[j],
            0,
            coord[center[i]][2] - coord[center[j]][2],
            exp_array[i],
            exp_array[j],
        )
        # ! derecha
        pz = 2.0 * exp_array[j] * phi.E(
            lz[i],
            lz[j] + 1,
            0,
            coord[center[i]][2] - coord[center[j]][2],
            exp_array[i],
            exp_array[j],
        ) - lz[j] * phi.E(
            lz[i],
            lz[j] - 1,
            0,
            coord[center[i]][2] - coord[center[j]][2],
            exp_array[i],
            exp_array[j],
        )

        # ! Terms of dxxLx
        xxd = (
            4.0
            * exp_array[i]
            * exp_array[i]
            * (
                phi.E(
                    lx[i] + 2,
                    lx[j],
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
            )
            - 2.0 * exp_array[i] * (2.0 * lx[i] + 1.0) * e0ij
            + lx[i]
            * (lx[i] - 1.0)
            * (
                phi.E(
                    lx[i] - 2,
                    lx[j],
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )

        dxxlx = xxd * ((e1kl + Ypg * e0kl) * pz - (e1mn + Zpg * e0mn) * py)

        # ! Terms of dyyLx
        ########### Dipole
        # E_1^kl + Ypc*E_0^kl
        d_skl_1 = 2.0 * exp_array[i] * (2.0 * ly[i] + 1) * (e1kl + Ypg * e0kl)
        # E_1^k-2l + Ypc*E_0^k-2l
        d_skm2l_1 = (
            ly[i]
            * (ly[i] - 1)
            * (
                phi.E(
                    ly[i] - 2,
                    ly[j],
                    1,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
                + Ypg
                * phi.E(
                    ly[i] - 2,
                    ly[j],
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        # E_1^k+2l + Ypc*E_0^k+2l
        d_skt2l_1 = (
            4.0
            * exp_array[i]
            * exp_array[i]
            * (
                phi.E(
                    ly[i] + 2,
                    ly[j],
                    1,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
                + Ypg
                * phi.E(
                    ly[i] + 2,
                    ly[j],
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        ####################### Derivatives
        # E_0^kl+1 - l*E_0^kl-1
        d_dkl_1 = (
            2.0
            * exp_array[i]
            * (2.0 * ly[i] + 1.0)
            * (
                2.0
                * exp_array[j]
                * phi.E(
                    ly[i],
                    ly[j] + 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
                - ly[j]
                * phi.E(
                    ly[i],
                    ly[j] - 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        # E_0^k-2l+1 - l*E_0^k-2l-1
        d_dkm2l_1 = (
            ly[i]
            * (ly[i] - 1.0)
            * (
                2.0
                * exp_array[j]
                * phi.E(
                    ly[i] - 2,
                    ly[j] + 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
                - ly[j]
                * phi.E(
                    ly[i] - 2,
                    ly[j] - 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        # E_0^k+2l+1 - l*E_0^k+2l-1
        d_dkt2l_1 = (
            4.0
            * exp_array[i]
            * exp_array[i]
            * (
                2.0
                * exp_array[j]
                * phi.E(
                    ly[i] + 2,
                    ly[j] + 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
                - ly[j]
                * phi.E(
                    ly[i] + 2,
                    ly[j] - 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )

        dyylx = (
            pz * (d_skm2l_1 - d_skl_1 + d_skt2l_1)
            + (e1mn + Zpg * e0mn) * (d_dkl_1 - d_dkm2l_1 - d_dkt2l_1)
        ) * e0ij

        # ! Terms of dzzLx
        # E_1^mn + Zpc*E_0^mn
        d_smn_1 = 2.0 * exp_array[i] * (2.0 * lz[i] + 1) * (e1mn + Zpg * e0mn)
        # E_1^m-2n + Zpc*E_0^m-2n
        d_smm2n_1 = (
            lz[i]
            * (lz[i] - 1)
            * (
                phi.E(
                    lz[i] - 2,
                    lz[j],
                    1,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
                + Zpg
                * phi.E(
                    lz[i] - 2,
                    lz[j],
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        # E_1^m+2n + Zpc*E_0^m+2n
        d_smt2n_1 = (
            4.0
            * exp_array[i]
            * exp_array[i]
            * (
                phi.E(
                    lz[i] + 2,
                    lz[j],
                    1,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
                + Zpg
                * phi.E(
                    lz[i] + 2,
                    lz[j],
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        ####################### Derivatives
        # E_0^kl+1 - l*E_0^kl-1
        d_dmn_1 = (
            2.0
            * exp_array[i]
            * (2.0 * lz[i] + 1.0)
            * (
                2.0
                * exp_array[j]
                * phi.E(
                    lz[i],
                    lz[j] + 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
                - lz[j]
                * phi.E(
                    lz[i],
                    lz[j] - 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        # E_0^k-2l+1 - l*E_0^k-2l-1
        d_dmm2n_1 = (
            lz[i]
            * (lz[i] - 1.0)
            * (
                2.0
                * exp_array[j]
                * phi.E(
                    lz[i] - 2,
                    lz[j] + 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
                - lz[j]
                * phi.E(
                    lz[i] - 2,
                    lz[j] - 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        # E_0^k+2l+1 - l*E_0^k+2l-1
        d_dmt2n_1 = (
            4.0
            * exp_array[i]
            * exp_array[i]
            * (
                2.0
                * exp_array[j]
                * phi.E(
                    lz[i] + 2,
                    lz[j] + 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
                - lz[j]
                * phi.E(
                    lz[i] + 2,
                    lz[j] - 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )

        dzzlx = (
            py * (d_smn_1 - d_smm2n_1 - d_smt2n_1)
            + (e1kl + Ypg * e0kl) * (d_dmm2n_1 - d_dmn_1 + d_dmt2n_1)
        ) * e0ij

        # * nabla Real{Lx} + Real{Lx} nabla

        intLx[i, j] = (
            Norm[n[i]](exp_array[i])
            * Norm[n[j]](exp_array[j])
            * 0.25
            * (dxxlx + dyylx + dzzlx)
            * np.power(np.pi / (exp_array[i] + exp_array[j]), 1.5)
        )
        intLx[j, i] = -1.0 * intLx[i, j]

        ###########################
        # ! Terms of dyyLy
        yyd = (
            4.0
            * exp_array[i]
            * exp_array[i]
            * (
                phi.E(
                    ly[i] + 2,
                    ly[j],
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
            )
            - 2.0 * exp_array[i] * (2.0 * ly[i] + 1.0) * e0kl
            + ly[i]
            * (ly[i] - 1.0)
            * (
                phi.E(
                    ly[i] - 2,
                    ly[j],
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )

        dyyly = yyd * ((e1mn + Zpg * e0mn) * px - (e1ij + Xpg * e0ij) * pz)

        # ! Terms of dzzLy
        ########### Dipole
        # E_1^kl + Ypc*E_0^kl
        d_skl_1 = 2.0 * exp_array[i] * (2.0 * lz[i] + 1) * (e1mn + Zpg * e0mn)
        # E_1^k-2l + Ypc*E_0^k-2l
        d_skm2l_1 = (
            lz[i]
            * (lz[i] - 1)
            * (
                phi.E(
                    lz[i] - 2,
                    lz[j],
                    1,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
                + Zpg
                * phi.E(
                    lz[i] - 2,
                    lz[j],
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        # E_1^k+2l + Ypc*E_0^k+2l
        d_skt2l_1 = (
            4.0
            * exp_array[i]
            * exp_array[i]
            * (
                phi.E(
                    lz[i] + 2,
                    lz[j],
                    1,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
                + Zpg
                * phi.E(
                    lz[i] + 2,
                    lz[j],
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        ####################### Derivatives
        # E_0^kl+1 - l*E_0^kl-1
        d_dkl_1 = (
            2.0
            * exp_array[i]
            * (2.0 * lz[i] + 1.0)
            * (
                2.0
                * exp_array[j]
                * phi.E(
                    lz[i],
                    lz[j] + 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
                - lz[j]
                * phi.E(
                    lz[i],
                    lz[j] - 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        # E_0^k-2l+1 - l*E_0^k-2l-1
        d_dkm2l_1 = (
            lz[i]
            * (lz[i] - 1.0)
            * (
                2.0
                * exp_array[j]
                * phi.E(
                    lz[i] - 2,
                    lz[j] + 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
                - lz[j]
                * phi.E(
                    lz[i] - 2,
                    lz[j] - 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        # E_0^k+2l+1 - l*E_0^k+2l-1
        d_dkt2l_1 = (
            4.0
            * exp_array[i]
            * exp_array[i]
            * (
                2.0
                * exp_array[j]
                * phi.E(
                    lz[i] + 2,
                    lz[j] + 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
                - lz[j]
                * phi.E(
                    lz[i] + 2,
                    lz[j] - 1,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )

        dzzly = (
            px * (d_skm2l_1 - d_skl_1 + d_skt2l_1)
            + (e1ij + Xpg * e0ij) * (d_dkl_1 - d_dkm2l_1 - d_dkt2l_1)
        ) * e0kl
        # ! Terms of dxxLy
        # E_1^mn + Zpc*E_0^mn
        d_smn_1 = 2.0 * exp_array[i] * (2.0 * lx[i] + 1) * (e1ij + Xpg * e0ij)
        # E_1^m-2n + Zpc*E_0^m-2n
        d_smm2n_1 = (
            lx[i]
            * (lx[i] - 1)
            * (
                phi.E(
                    lx[i] - 2,
                    lx[j],
                    1,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
                + Xpg
                * phi.E(
                    lx[i] - 2,
                    lx[j],
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        # E_1^m+2n + Zpc*E_0^m+2n
        d_smt2n_1 = (
            4.0
            * exp_array[i]
            * exp_array[i]
            * (
                phi.E(
                    lx[i] + 2,
                    lx[j],
                    1,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
                + Xpg
                * phi.E(
                    lx[i] + 2,
                    lx[j],
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        ####################### Derivatives
        # E_0^kl+1 - l*E_0^kl-1
        d_dmn_1 = (
            2.0
            * exp_array[i]
            * (2.0 * lx[i] + 1.0)
            * (
                2.0
                * exp_array[j]
                * phi.E(
                    lx[i],
                    lx[j] + 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
                - lx[j]
                * phi.E(
                    lx[i],
                    lx[j] - 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        # E_0^k-2l+1 - l*E_0^k-2l-1
        d_dmm2n_1 = (
            lx[i]
            * (lx[i] - 1.0)
            * (
                2.0
                * exp_array[j]
                * phi.E(
                    lx[i] - 2,
                    lx[j] + 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
                - lx[j]
                * phi.E(
                    lx[i] - 2,
                    lx[j] - 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        # E_0^k+2l+1 - l*E_0^k+2l-1
        d_dmt2n_1 = (
            4.0
            * exp_array[i]
            * exp_array[i]
            * (
                2.0
                * exp_array[j]
                * phi.E(
                    lx[i] + 2,
                    lx[j] + 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
                - lx[j]
                * phi.E(
                    lx[i] + 2,
                    lx[j] - 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )

        dxxly = (
            pz * (d_smn_1 - d_smm2n_1 - d_smt2n_1)
            + (e1mn + Zpg * e0mn) * (d_dmm2n_1 - d_dmn_1 + d_dmt2n_1)
        ) * e0kl
        # * nabla Real{Ly} + Real{Ly} nabla

        intLy[i, j] = (
            Norm[n[i]](exp_array[i])
            * Norm[n[j]](exp_array[j])
            * 0.25
            * (dxxly + dyyly + dzzly)
            * np.power(np.pi / (exp_array[i] + exp_array[j]), 1.5)
        )
        intLy[j, i] = -1.0 * intLy[i, j]

        ###########################
        # ! Terms of dzzLz
        zzd = (
            4.0
            * exp_array[i]
            * exp_array[i]
            * (
                phi.E(
                    lz[i] + 2,
                    lz[j],
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
            )
            - 2.0 * exp_array[i] * (2.0 * lz[i] + 1.0) * e0mn
            + lz[i]
            * (lz[i] - 1.0)
            * (
                phi.E(
                    lz[i] - 2,
                    lz[j],
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )

        dzzlz = zzd * ((e1ij + Xpg * e0ij) * py - (e1kl + Ypg * e0kl) * px)

        # ! Terms of dxxLz
        ########### Dipole

        d_skl_1 = 2.0 * exp_array[i] * (2.0 * lx[i] + 1) * (e1ij + Xpg * e0ij)

        d_skm2l_1 = (
            lx[i]
            * (lx[i] - 1)
            * (
                phi.E(
                    lx[i] - 2,
                    lx[j],
                    1,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
                + Xpg
                * phi.E(
                    lx[i] - 2,
                    lx[j],
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        # E_1^k+2l + Ypc*E_0^k+2l
        d_skt2l_1 = (
            4.0
            * exp_array[i]
            * exp_array[i]
            * (
                phi.E(
                    lx[i] + 2,
                    lx[j],
                    1,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
                + Xpg
                * phi.E(
                    lx[i] + 2,
                    lx[j],
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        ####################### Derivatives
        # E_0^kl+1 - l*E_0^kl-1
        d_dkl_1 = (
            2.0
            * exp_array[i]
            * (2.0 * lx[i] + 1.0)
            * (
                2.0
                * exp_array[j]
                * phi.E(
                    lx[i],
                    lx[j] + 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
                - lx[j]
                * phi.E(
                    lx[i],
                    lx[j] - 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )

        d_dkm2l_1 = (
            lx[i]
            * (lx[i] - 1.0)
            * (
                2.0
                * exp_array[j]
                * phi.E(
                    lx[i] - 2,
                    lx[j] + 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
                - lx[j]
                * phi.E(
                    lx[i] - 2,
                    lx[j] - 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )

        d_dkt2l_1 = (
            4.0
            * exp_array[i]
            * exp_array[i]
            * (
                2.0
                * exp_array[j]
                * phi.E(
                    lx[i] + 2,
                    lx[j] + 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
                - lx[j]
                * phi.E(
                    lx[i] + 2,
                    lx[j] - 1,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )

        dxxlz = (
            py * (d_skm2l_1 - d_skl_1 + d_skt2l_1)
            + (e1kl + Ypg * e0kl) * (d_dkl_1 - d_dkm2l_1 - d_dkt2l_1)
        ) * e0mn

        # ! Terms of dyyLz

        d_smn_1 = 2.0 * exp_array[i] * (2.0 * ly[i] + 1) * (e1kl + Ypg * e0kl)

        d_smm2n_1 = (
            ly[i]
            * (ly[i] - 1)
            * (
                phi.E(
                    ly[i] - 2,
                    ly[j],
                    1,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
                + Ypg
                * phi.E(
                    ly[i] - 2,
                    ly[j],
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )

        d_smt2n_1 = (
            4.0
            * exp_array[i]
            * exp_array[i]
            * (
                phi.E(
                    ly[i] + 2,
                    ly[j],
                    1,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
                + Ypg
                * phi.E(
                    ly[i] + 2,
                    ly[j],
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )
        ####################### Derivatives

        d_dmn_1 = (
            2.0
            * exp_array[i]
            * (2.0 * ly[i] + 1.0)
            * (
                2.0
                * exp_array[j]
                * phi.E(
                    ly[i],
                    ly[j] + 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
                - ly[j]
                * phi.E(
                    ly[i],
                    ly[j] - 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )

        d_dmm2n_1 = (
            ly[i]
            * (ly[i] - 1.0)
            * (
                2.0
                * exp_array[j]
                * phi.E(
                    ly[i] - 2,
                    ly[j] + 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
                - ly[j]
                * phi.E(
                    ly[i] - 2,
                    ly[j] - 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )

        d_dmt2n_1 = (
            4.0
            * exp_array[i]
            * exp_array[i]
            * (
                2.0
                * exp_array[j]
                * phi.E(
                    ly[i] + 2,
                    ly[j] + 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
                - ly[j]
                * phi.E(
                    ly[i] + 2,
                    ly[j] - 1,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
            )
        )

        dyylz = (
            px * (d_smn_1 - d_smm2n_1 - d_smt2n_1)
            + (e1ij + Xpg * e0ij) * (d_dmm2n_1 - d_dmn_1 + d_dmt2n_1)
        ) * e0mn
        # * nabla Real{Lz} + Real{Lz} nabla

        intLz[i, j] = (
            Norm[n[i]](exp_array[i])
            * Norm[n[j]](exp_array[j])
            * 0.25
            * (dxxlz + dyylz + dzzlz)
            * np.power(np.pi / (exp_array[i] + exp_array[j]), 1.5)
        )
        intLz[j, i] = -1.0 * intLz[i, j]

        if (
            output > 10
            and np.abs(intLz[j, i]) > 1e-10  # and j == 7 and i == 0
        ):  # :
            print(
                "int [",
                j + 1,
                ",",
                i + 1,
                "] : ",
                intLz[j, i],
            )

print(" time [s]: ", -start + time.time())
