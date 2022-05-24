from numpy import exp
import numpy as np
import time
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

intLx = np.zeros((2, total_nprim, total_nprim), dtype=float)
intLy = np.zeros((2, total_nprim, total_nprim), dtype=float)
intLz = np.zeros((2, total_nprim, total_nprim), dtype=float)
output = 9

# !Note: Rgaugeo is put on in the origen
for k in range(2):
    if output > 10:
        print(" **** Atom ", k + 1, " **** \n")
    for i in range(total_nprim):

        for j in range(total_nprim):

            ydz = 2.0 * exp[j] * nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j],
                ly[j],
                lz[j] + 1,
                0,
                1,
                0,
                exp[i],
                exp[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                coord[k][0],
                coord[k][1],
                coord[k][2],
            ) - lz[j] * nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j],
                ly[j],
                lz[j] - 1,
                0,
                1,
                0,
                exp[i],
                exp[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                coord[k][0],
                coord[k][1],
                coord[k][2],
            )

            zdy = ( 
                2.0 * exp[j] * nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j],
                ly[j] + 1,
                lz[j],
                0,
                0,
                1,
                exp[i],
                exp[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                coord[k][0],
                coord[k][1],
                coord[k][2],
            ) - 
                ly[j] * nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j],
                ly[j] - 1,
                lz[j],
                0,
                0,
                1,
                exp[i],
                exp[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                coord[k][0],
                coord[k][1],
                coord[k][2],
            ))

            xdz = 2.0 * exp[j] * nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j],
                ly[j],
                lz[j] + 1,
                1,
                0,
                0,
                exp[i],
                exp[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                coord[k][0],
                coord[k][1],
                coord[k][2],
            ) - lz[j] * nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j],
                ly[j],
                lz[j] - 1,
                1,
                0,
                0,
                exp[i],
                exp[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                coord[k][0],
                coord[k][1],
                coord[k][2],
            )

            zdx = 2.0 * exp[j] * nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j] + 1,
                ly[j],
                lz[j],
                0,
                0,
                1,
                exp[i],
                exp[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                coord[k][0],
                coord[k][1],
                coord[k][2],
            ) - lx[j] * nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j] - 1,
                ly[j],
                lz[j],
                0,
                0,
                1,
                exp[i],
                exp[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                coord[k][0],
                coord[k][1],
                coord[k][2],
            )

            xdy = 2.0 * exp[j] * nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j],
                ly[j] + 1,
                lz[j],
                1,
                0,
                0,
                exp[i],
                exp[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                coord[k][0],
                coord[k][1],
                coord[k][2],
            ) - ly[j] * nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j],
                ly[j] - 1,
                lz[j],
                1,
                0,
                0,
                exp[i],
                exp[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                coord[k][0],
                coord[k][1],
                coord[k][2],
            )

            ydx = 2.0 * exp[j] * nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j] + 1,
                ly[j],
                lz[j],
                0,
                1,
                0,
                exp[i],
                exp[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                coord[k][0],
                coord[k][1],
                coord[k][2],
            ) - lx[j] * nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j] - 1,
                ly[j],
                lz[j],
                0,
                1,
                0,
                exp[i],
                exp[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                coord[k][0],
                coord[k][1],
                coord[k][2],
            )

            # ! Terms of dxxLx/rk^3
            xxdydz = (
                4.0
                * exp[i]
                * exp[i]
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lz[j]
                    * nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] - 1,
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
                - 2.0 * exp[i] * (2.0 * lx[i] + 1.0) * ydz  # e0ij
                + lx[i]
                * (lx[i] - 1.0)
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lz[j]
                    * nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] - 1,
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )
            # print("dxxydz",xxdydz)

            xxdzdy = (
                4.0
                * exp[i]
                * exp[i]
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - ly[j]
                    * nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] - 1,
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            - 2.0 * exp[i] * (2.0 * lx[i] + 1.0) * zdy
            + lx[i]
                * (lx[i] - 1.0)
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - ly[j]
                    * nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] - 1,
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )
            # print("dxxzdy",xxdzdy)

            dxxlx = xxdydz - xxdzdy

            # ! Terms of dyyLx
            ########### Dipole

            d_skl_1 = (
                2.0 * exp[i] * (2.0 * ly[i] + 1) * ydz
            )  # (e1kl + Ypg * e0kl)

            d_skm2l_1 = (
                ly[i]
                * (ly[i] - 1)
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lz[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] - 1,
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )
            # E_1^k+2l + Ypc*E_0^k+2l
            d_skt2l_1 = (
                4.0
                * exp[i]
                * exp[i]
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lz[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] - 1,
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
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
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - ly[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] - 1,
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
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
                    * nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - ly[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j],
                        ly[j] - 1,
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
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
                    * nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - ly[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j],
                        ly[j] - 1,
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )

            dyylx = (d_skm2l_1 - d_skl_1 + d_skt2l_1) - (
                d_dkm2l_1 - d_dkl_1 + d_dkt2l_1
            )

            # ! Terms of dzzLx
            # E_1^mn + Zpc*E_0^mn
            d_smn_1 = (
                2.0
                * exp[i]
                * (2.0 * lz[i] + 1)
                * ydz  # (e1mn + Zpg * e0mn)
            )
            # E_1^m-2n + Zpc*E_0^m-2n
            d_smm2n_1 = (
                lz[i]
                * (lz[i] - 1)
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lz[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j],
                        ly[j],
                        lz[j] - 1,
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )
            # E_1^m+2n + Zpc*E_0^m+2n
            d_smt2n_1 = (
                4.0
                * exp[i]
                * exp[i]
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lz[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j],
                        ly[j],
                        lz[j] - 1,
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
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
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - ly[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] - 1,
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
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
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - ly[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j],
                        ly[j] - 1,
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
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
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - ly[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j],
                        ly[j] - 1,
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )

            dzzlx = (d_smm2n_1 - d_smn_1 + d_smt2n_1) - (
                d_dmm2n_1 - d_dmn_1 + d_dmt2n_1
            )

            # * nabla Real{Lx} + Real{Lx} nabla
            #print("ydz-zdy",i,j,dxxlx + dyylx + dzzlx)

            intLx[k, j, i] = (
                -Norm[lx[i] + ly[i] + lz[i]](exp[i])
                * Norm[lx[j] + ly[j] + lz[j]](exp[j])
                * (dxxlx + dyylx + dzzlx)
                * 0.5
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
            )

            # ! Terms of dyyLy/rk^3
            yydzdx = (
                4.0
                * exp[i]
                * exp[i]
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lx[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j] - 1,
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
                - 2.0 * exp[i] * (2.0 * ly[i] + 1.0) * zdx  # e0ij
                + ly[i]
                * (ly[i] - 1.0)
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lx[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j] - 1,
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )

            yydxdz = (
                4.0
                * exp[i]
                * exp[i]
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lz[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] - 1,
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            ) - 2.0 * exp[i] * (2.0 * ly[i] + 1.0) * xdz
            +(
                ly[i]
                * (ly[i] - 1.0)
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lz[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] - 1,
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )

            dyyly = yydzdx - yydxdz

            # ! Terms of dxxLy
            ########### Dipole

            d_skl_1 = (
                2.0 * exp[i] * (2.0 * lx[i] + 1) * zdx
            )  # (e1kl + Ypg * e0kl)

            d_skm2l_1 = (
                lx[i]
                * (lx[i] - 1)
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lx[j]
                    * nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j] - 1,
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )
            # E_1^k+2l + Ypc*E_0^k+2l
            d_skt2l_1 = (
                4.0
                * exp[i]
                * exp[i]
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lx[j]
                    * nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j] - 1,
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
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
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lz[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] - 1,
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )
            # E_0^k-2l+1 - l*E_0^k-2l-1
            d_dkm2l_1 = (
                lx[i]
                * (lx[i] - 1.0)
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lz[j]
                    * nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] - 1,
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
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
                    * nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lz[j]
                    * nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] - 1,
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )

            dxxly = (d_skm2l_1 - d_skl_1 + d_skt2l_1) - (
                d_dkm2l_1 - d_dkl_1 + d_dkt2l_1
            )

            # ! Terms of dzzLy
            # E_1^mn + Zpc*E_0^mn
            d_smn_1 = 2.0 * exp[i] * (2.0 * lz[i] + 1) * zdx
            # E_1^m-2n + Zpc*E_0^m-2n
            d_smm2n_1 = (
                lz[i]
                * (lz[i] - 1)
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lx[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j] - 1,
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )
            # E_1^m+2n + Zpc*E_0^m+2n
            d_smt2n_1 = (
                4.0
                * exp[i]
                * exp[i]
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lx[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j] - 1,
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
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
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lz[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] - 1,
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
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
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lz[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j],
                        ly[j],
                        lz[j] - 1,
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
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
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lz[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j],
                        ly[j],
                        lz[j] - 1,
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )

            dzzly = (d_smm2n_1 - d_smn_1 + d_smt2n_1) - (
                d_dmm2n_1 - d_dmn_1 + d_dmt2n_1
            )

            # * nabla Real{Ly} + Real{Ly} nabla

            intLy[k, j, i] = (
                -Norm[lx[i] + ly[i] + lz[i]](exp[i])
                * Norm[lx[j] + ly[j] + lz[j]](exp[j])
                * (dxxly + dyyly + dzzly)
                * 0.5
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
            )

            # ! Terms of dyyLy/rk^3
            zzdxdy = (
                4.0
                * exp[i]
                * exp[i]
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - ly[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j],
                        ly[j] - 1,
                        lz[j],
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
                - 2.0 * exp[i] * (2.0 * lz[i] + 1.0) * xdy
                + lz[i]
                * (lz[i] - 1.0)
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - ly[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j],
                        ly[j] - 1,
                        lz[j],
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )

            zzdydx = (
                4.0
                * exp[i]
                * exp[i]
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lx[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j] - 1,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            ) - 2.0 * exp[i] * (2.0 * lz[i] + 1.0) * ydx
            +(
                lz[i]
                * (lz[i] - 1.0)
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lx[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j] - 1,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )

            dzzlz = zzdxdy - zzdydx

            # ! Terms of dxxLy
            ########### Dipole

            d_skl_1 = (
                2.0 * exp[i] * (2.0 * lx[i] + 1) * xdy
            )  # (e1kl + Ypg * e0kl)

            d_skm2l_1 = (
                lx[i]
                * (lx[i] - 1)
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - ly[j]
                    * nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] - 1,
                        lz[j],
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )
            # E_1^k+2l + Ypc*E_0^k+2l
            d_skt2l_1 = (
                4.0
                * exp[i]
                * exp[i]
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - ly[j]
                    * nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] - 1,
                        lz[j],
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
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
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lx[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] - 1,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )
            # E_0^k-2l+1 - l*E_0^k-2l-1
            d_dkm2l_1 = (
                lx[i]
                * (lx[i] - 1.0)
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lx[j]
                    * nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j] - 1,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
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
                    * nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lx[j]
                    * nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j] - 1,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )

            dxxlz = (d_skm2l_1 - d_skl_1 + d_skt2l_1) - (
                d_dkm2l_1 - d_dkl_1 + d_dkt2l_1
            )

            # ! Terms of dyyLz
            # E_1^mn + Zpc*E_0^mn
            d_smn_1 = 2.0 * exp[i] * (2.0 * ly[i] + 1) * xdy

            # E_1^m-2n + Zpc*E_0^m-2n
            d_smm2n_1 = (
                ly[i]
                * (ly[i] - 1)
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - ly[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j],
                        ly[j] - 1,
                        lz[j],
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )

            # E_1^m+2n + Zpc*E_0^m+2n
            d_smt2n_1 = (
                4.0
                * exp[i]
                * exp[i]
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - ly[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j],
                        ly[j] - 1,
                        lz[j],
                        1,
                        0,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )

            ####################### Derivatives
            # E_0^kl+1 - l*E_0^kl-1
            d_dmn_1 = (
                2.0
                * exp[i]
                * (2.0 * ly[i] + 1.0)
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lx[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] - 1,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )
            # E_0^k-2l+1 - l*E_0^k-2l-1
            d_dmm2n_1 = (
                ly[i]
                * (ly[i] - 1.0)
                * (
                    2.0
                    * exp[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lx[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j] - 1,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
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
                    * nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                    - lx[j]
                    * nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j] - 1,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        coord[k][0],
                        coord[k][1],
                        coord[k][2],
                    )
                )
            )

            dyylz = (d_smm2n_1 - d_smn_1 + d_smt2n_1) - (
                d_dmm2n_1 - d_dmn_1 + d_dmt2n_1
            )
            # * nabla Real{Ly} + Real{Ly} nabla

            intLz[k, j, i] = (
                -Norm[lx[i] + ly[i] + lz[i]](exp[i])
                * Norm[lx[j] + ly[j] + lz[j]](exp[j])
                * (dxxlz + dyylz + dzzlz)
                * 0.5
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
            )

            if output > 1 and abs(intLx[k, i, j]) > 0.002:
                print(
                    "int [",
                    i + 1,
                    ",",
                    j + 1,
                    "] : ",
                    intLx[k, i, j],
                )


print(" time [s]: ", -start + time())
