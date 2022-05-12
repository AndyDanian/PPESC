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
intMASS = np.zeros((total_nprim, total_nprim), dtype=float)
output = 11

c = 137.0359998
a2 = 1.0 / (c * c)

#! Note: It's neccesary to multiplicate by 5.0/3.0 to get DALTON's values
const_mass = a2 / 8.0 * 5.0 / 3.0
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

        sij = phi.E(
            lx[i],
            lx[j],
            0,
            coord[center[i]][0] - coord[center[j]][0],
            exp_array[i],
            exp_array[j],
        )

        skl = phi.E(
            ly[i],
            ly[j],
            0,
            coord[center[i]][1] - coord[center[j]][1],
            exp_array[i],
            exp_array[j],
        )

        smn = phi.E(
            lz[i],
            lz[j],
            0,
            coord[center[i]][2] - coord[center[j]][2],
            exp_array[i],
            exp_array[j],
        )

        # dxxdxx
        dxxsidxxsj = (
            16.0
            * exp_array[i]
            * exp_array[i]
            * exp_array[j]
            * exp_array[j]
            * (
                phi.E(
                    lx[i] + 2,
                    lx[j] + 2,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
            )
            + 4.0
            * exp_array[i]
            * exp_array[j]
            * (2.0 * lx[i] + 1.0)
            * (2.0 * lx[j] + 1.0)
            * sij
            + lx[i]
            * (lx[i] - 1.0)
            * lx[j]
            * (lx[j] - 1.0)
            * (
                phi.E(
                    lx[i] - 2,
                    lx[j] - 2,
                    0,
                    coord[center[i]][0] - coord[center[j]][0],
                    exp_array[i],
                    exp_array[j],
                )
            )
            + 4.0
            * exp_array[i]
            * exp_array[i]
            * lx[j]
            * (lx[j] - 1.0)
            * phi.E(
                lx[i] + 2,
                lx[j] - 2,
                0,
                coord[center[i]][0] - coord[center[j]][0],
                exp_array[i],
                exp_array[j],
            )
            + 4.0
            * exp_array[j]
            * exp_array[j]
            * lx[i]
            * (lx[i] - 1.0)
            * phi.E(
                lx[i] - 2,
                lx[j] + 2,
                0,
                coord[center[i]][0] - coord[center[j]][0],
                exp_array[i],
                exp_array[j],
            )
            - 8.0
            * exp_array[i]
            * exp_array[i]
            * exp_array[j]
            * (2.0 * lx[j] + 1.0)
            * phi.E(
                lx[i] + 2,
                lx[j],
                0,
                coord[center[i]][0] - coord[center[j]][0],
                exp_array[i],
                exp_array[j],
            )
            - 8.0
            * exp_array[j]
            * exp_array[j]
            * exp_array[i]
            * (2.0 * lx[i] + 1.0)
            * phi.E(
                lx[i],
                lx[j] + 2,
                0,
                coord[center[i]][0] - coord[center[j]][0],
                exp_array[i],
                exp_array[j],
            )
            - 2.0
            * exp_array[i]
            * (2.0 * lx[i] + 1)
            * lx[j]
            * (lx[j] - 1)
            * phi.E(
                lx[i],
                lx[j] - 2,
                0,
                coord[center[i]][0] - coord[center[j]][0],
                exp_array[i],
                exp_array[j],
            )
            - 2.0
            * exp_array[j]
            * (2.0 * lx[j] + 1)
            * lx[i]
            * (lx[i] - 1)
            * phi.E(
                lx[i] - 2,
                lx[j],
                0,
                coord[center[i]][0] - coord[center[j]][0],
                exp_array[i],
                exp_array[j],
            )
        )

        # dyydyy,
        dyyskdyysl = (
            16.0
            * exp_array[i]
            * exp_array[i]
            * exp_array[j]
            * exp_array[j]
            * (
                phi.E(
                    ly[i] + 2,
                    ly[j] + 2,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
            )
            + 4.0
            * exp_array[i]
            * exp_array[j]
            * (2.0 * ly[i] + 1.0)
            * (2.0 * ly[j] + 1.0)
            * skl
            + ly[i]
            * (ly[i] - 1.0)
            * ly[j]
            * (ly[j] - 1.0)
            * (
                phi.E(
                    ly[i] - 2,
                    ly[j] - 2,
                    0,
                    coord[center[i]][1] - coord[center[j]][1],
                    exp_array[i],
                    exp_array[j],
                )
            )
            + 4.0
            * exp_array[i]
            * exp_array[i]
            * ly[j]
            * (ly[j] - 1.0)
            * phi.E(
                ly[i] + 2,
                ly[j] - 2,
                0,
                coord[center[i]][1] - coord[center[j]][1],
                exp_array[i],
                exp_array[j],
            )
            + 4.0
            * exp_array[j]
            * exp_array[j]
            * ly[i]
            * (ly[i] - 1.0)
            * phi.E(
                ly[i] - 2,
                ly[j] + 2,
                0,
                coord[center[i]][1] - coord[center[j]][1],
                exp_array[i],
                exp_array[j],
            )
            - 8.0
            * exp_array[i]
            * exp_array[i]
            * exp_array[j]
            * (2.0 * ly[j] + 1.0)
            * phi.E(
                ly[i] + 2,
                ly[j],
                0,
                coord[center[i]][1] - coord[center[j]][1],
                exp_array[i],
                exp_array[j],
            )
            - 8.0
            * exp_array[j]
            * exp_array[j]
            * exp_array[i]
            * (2.0 * ly[i] + 1.0)
            * phi.E(
                ly[i],
                ly[j] + 2,
                0,
                coord[center[i]][1] - coord[center[j]][1],
                exp_array[i],
                exp_array[j],
            )
            - 2.0
            * exp_array[i]
            * (2.0 * ly[i] + 1)
            * ly[j]
            * (ly[j] - 1)
            * phi.E(
                ly[i],
                ly[j] - 2,
                0,
                coord[center[i]][1] - coord[center[j]][1],
                exp_array[i],
                exp_array[j],
            )
            - 2.0
            * exp_array[j]
            * (2.0 * ly[j] + 1)
            * ly[i]
            * (ly[i] - 1)
            * phi.E(
                ly[i] - 2,
                ly[j],
                0,
                coord[center[i]][1] - coord[center[j]][1],
                exp_array[i],
                exp_array[j],
            )
        )

        # dzzdzz,
        dzzsmdzzsn = (
            16.0
            * exp_array[i]
            * exp_array[i]
            * exp_array[j]
            * exp_array[j]
            * (
                phi.E(
                    lz[i] + 2,
                    lz[j] + 2,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
            )
            + 4.0
            * exp_array[i]
            * exp_array[j]
            * (2.0 * lz[i] + 1.0)
            * (2.0 * lz[j] + 1.0)
            * smn
            + lz[i]
            * (lz[i] - 1.0)
            * lz[j]
            * (lz[j] - 1.0)
            * (
                phi.E(
                    lz[i] - 2,
                    lz[j] - 2,
                    0,
                    coord[center[i]][2] - coord[center[j]][2],
                    exp_array[i],
                    exp_array[j],
                )
            )
            + 4.0
            * exp_array[i]
            * exp_array[i]
            * lz[j]
            * (lz[j] - 1.0)
            * phi.E(
                lz[i] + 2,
                lz[j] - 2,
                0,
                coord[center[i]][2] - coord[center[j]][2],
                exp_array[i],
                exp_array[j],
            )
            + 4.0
            * exp_array[j]
            * exp_array[j]
            * lz[i]
            * (lz[i] - 1.0)
            * phi.E(
                lz[i] - 2,
                lz[j] + 2,
                0,
                coord[center[i]][2] - coord[center[j]][2],
                exp_array[i],
                exp_array[j],
            )
            - 8.0
            * exp_array[i]
            * exp_array[i]
            * exp_array[j]
            * (2.0 * lz[j] + 1.0)
            * phi.E(
                lz[i] + 2,
                lz[j],
                0,
                coord[center[i]][2] - coord[center[j]][2],
                exp_array[i],
                exp_array[j],
            )
            - 8.0
            * exp_array[j]
            * exp_array[j]
            * exp_array[i]
            * (2.0 * lz[i] + 1.0)
            * phi.E(
                lz[i],
                lz[j] + 2,
                0,
                coord[center[i]][2] - coord[center[j]][2],
                exp_array[i],
                exp_array[j],
            )
            - 2.0
            * exp_array[i]
            * (2.0 * lz[i] + 1)
            * lz[j]
            * (lz[j] - 1)
            * phi.E(
                lz[i],
                lz[j] - 2,
                0,
                coord[center[i]][2] - coord[center[j]][2],
                exp_array[i],
                exp_array[j],
            )
            - 2.0
            * exp_array[j]
            * (2.0 * lz[j] + 1)
            * lz[i]
            * (lz[i] - 1)
            * phi.E(
                lz[i] - 2,
                lz[j],
                0,
                coord[center[i]][2] - coord[center[j]][2],
                exp_array[i],
                exp_array[j],
            )
        )

        intMASS[i, j] = intMASS[j, i] = (
            -Norm[n[i]](exp_array[i])
            * Norm[n[j]](exp_array[j])
            * const_mass
            * (
                dxxsidxxsj * skl * smn
                + sij * dyyskdyysl * smn
                + sij * skl * dzzsmdzzsn
            )
            * np.power(np.pi / (exp_array[i] + exp_array[j]), 1.5)
        )

        if output > 10:
            print("int [", i + 1, ",", j + 1, "] : ", intMASS[i, j])

print(" time [s]: ", -start + time.time())
