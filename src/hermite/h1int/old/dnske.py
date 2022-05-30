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
intNEFx = np.zeros((2, total_nprim, total_nprim), dtype=float)
intNEFy = np.zeros((2, total_nprim, total_nprim), dtype=float)
intNEFz = np.zeros((2, total_nprim, total_nprim), dtype=float)
output = 11

Rg = [0.0, 0.0, 1.4045523587]
for k in range(2):
    print("\n   ****Atom  ", k + 1, " ****\n")
    for i in range(total_nprim):

        for j in range(i, total_nprim):

            # Nuclear Diamagnetic Shielding
            # <phi|[nabla^2, (y-yg)(y-yk)/rk^3 + (z-zg)(z-zk)/rk^3]_+|phi> --> Eqs 9.932, 9.9.18-20

            A2x = (
                phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j] + 1,
                    lz[j],
                    0,
                    1,
                    0,
                    exp_array[i],
                    exp_array[j],
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
                + (coord[center[j]][1] - Rg[1])
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j],
                    lz[j],
                    0,
                    1,
                    0,
                    exp_array[i],
                    exp_array[j],
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
                + phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j],
                    lz[j] + 1,
                    0,
                    0,
                    1,
                    exp_array[i],
                    exp_array[j],
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
                + (coord[center[j]][2] - Rg[2])
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j],
                    lz[j],
                    0,
                    0,
                    1,
                    exp_array[i],
                    exp_array[j],
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

            A2y = (
                phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i],
                    lx[j] + 1,
                    ly[j],
                    lz[j],
                    1,
                    0,
                    0,
                    exp_array[i],
                    exp_array[j],
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
                + (coord[center[j]][0] - Rg[0])
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j],
                    lz[j],
                    1,
                    0,
                    0,
                    exp_array[i],
                    exp_array[j],
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
                + phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j],
                    lz[j] + 1,
                    0,
                    0,
                    1,
                    exp_array[i],
                    exp_array[j],
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
                + (coord[center[j]][2] - Rg[2])
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j],
                    lz[j],
                    0,
                    0,
                    1,
                    exp_array[i],
                    exp_array[j],
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

            A2z = (
                phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i],
                    lx[j] + 1,
                    ly[j],
                    lz[j],
                    1,
                    0,
                    0,
                    exp_array[i],
                    exp_array[j],
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
                + (coord[center[j]][0] - Rg[0])
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j],
                    lz[j],
                    1,
                    0,
                    0,
                    exp_array[i],
                    exp_array[j],
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
                + phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j] + 1,
                    lz[j],
                    0,
                    1,
                    0,
                    exp_array[i],
                    exp_array[j],
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
                + (coord[center[j]][1] - Rg[1])
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j],
                    lz[j],
                    0,
                    1,
                    0,
                    exp_array[i],
                    exp_array[j],
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

            # *** dx^2(y-yg)(y-yk)/r^3_k + (z-zg)(z-zk)/r^3_k ***
            dxxAx_a = (
                4.0
                * exp_array[i]
                * exp_array[i]
                * (
                    phi.nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            dxxAx_b = -2.0 * exp_array[i] * (2.0 * lx[i] + 1.0) * A2x

            dxxAx_c = (
                lx[i]
                * (lx[i] - 1.0)
                * (
                    phi.nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            # *** dx^2(y-yg)(y-yk)/r^3_k + (z-zg)(z-zk)/r^3_k ***
            dyyAx_a = (
                4.0
                * exp_array[i]
                * exp_array[i]
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            dyyAx_b = -2.0 * exp_array[i] * (2.0 * ly[i] + 1.0) * A2x

            dyyAx_c = (
                ly[i]
                * (ly[i] - 1.0)
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            # *** dz^2((y-yg)(y-yk)/r^3_k + (z-zg)(z-zk)/r^3_k) ***
            dzzAx_a = (
                4.0
                * exp_array[i]
                * exp_array[i]
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            dzzAx_b = -2.0 * exp_array[i] * (2.0 * lz[i] + 1.0) * A2x

            dzzAx_c = (
                lz[i]
                * (lz[i] - 1.0)
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            # *** ((y-yg)(y-yk)/r^3_k + (z-zg)(z-zk)/r^3_k)dx^2 ***
            Axdxx_a = (
                4.0
                * exp_array[j]
                * exp_array[j]
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 2,
                        ly[j] + 1,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 2,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 2,
                        ly[j],
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 2,
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            Axdxx_b = -2.0 * exp_array[j] * (2.0 * lx[j] + 1.0) * A2x

            Axdxx_c = (
                lx[j]
                * (lx[j] - 1.0)
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] - 2,
                        ly[j] + 1,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] - 2,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] - 2,
                        ly[j],
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] - 2,
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            # *** ((y-yg)(y-yk)/r^3_k + (z-zg)(z-zk)/r^3_k)dy^2 ***
            Axdyy_a = (
                4.0
                * exp_array[j]
                * exp_array[j]
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 3,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 2,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 2,
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 2,
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            Axdyy_b = -2.0 * exp_array[j] * (2.0 * ly[j] + 1.0) * A2x

            Axdyy_c = (
                ly[j]
                * (ly[j] - 1.0)
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] - 1,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] - 2,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] - 2,
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] - 2,
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            # *** ((y-yg)(y-yk)/r^3_k + (z-zg)(z-zk)/r^3_k)dz^2 ***
            Axdzz_a = (
                4.0
                * exp_array[j]
                * exp_array[j]
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j] + 2,
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 2,
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 3,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 2,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            Axdzz_b = -2.0 * exp_array[j] * (2.0 * lz[j] + 1.0) * A2x

            Axdzz_c = (
                lz[j]
                * (lz[j] - 1.0)
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j] - 2,
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] - 2,
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] - 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] - 2,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            # ! Falta poner el -1^e+g+f en la recurrencia de nuclear_attraction
            intNEFx[k, i, j] = intNEFx[k, j, i] = (
                -Norm[n[i]](exp_array[i])
                * Norm[n[j]](exp_array[j])
                * 2.0
                * np.pi
                / (exp_array[i] + exp_array[j])
                * (
                    dxxAx_a
                    + dxxAx_b
                    + dxxAx_c
                    + dyyAx_a
                    + dyyAx_b
                    + dyyAx_c
                    + dzzAx_a
                    + dzzAx_b
                    + dzzAx_c
                    + Axdxx_a
                    + Axdxx_b
                    + Axdxx_c
                    + Axdyy_a
                    + Axdyy_b
                    + Axdyy_c
                    + Axdzz_a
                    + Axdzz_b
                    + Axdzz_c
                )
                * 3.0
                / 4.0
            )

            # *** dxx((x-xg)(x-xk)/r^3_k + (y-yg)(y-yk)/r^3_k) ***
            dxxAy_a = (
                4.0
                * exp_array[i]
                * exp_array[i]
                * (
                    phi.nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            dxxAy_b = -2.0 * exp_array[i] * (2.0 * lx[i] + 1.0) * A2y

            dxxAy_c = (
                lx[i]
                * (lx[i] - 1.0)
                * (
                    phi.nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            # *** dyy((x-xg)(x-xk)/r^3_k + (y-yg)(y-yk)/r^3_k) ***
            dyyAy_a = (
                4.0
                * exp_array[i]
                * exp_array[i]
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            dyyAy_b = -2.0 * exp_array[i] * (2.0 * ly[i] + 1.0) * A2y

            dyyAy_c = (
                ly[i]
                * (ly[i] - 1.0)
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            # *** dzz((x-xg)(x-xk)/r^3_k + (y-yg)(y-yk)/r^3_k) ***
            dzzAy_a = (
                4.0
                * exp_array[i]
                * exp_array[i]
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j],
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            dzzAy_b = -2.0 * exp_array[i] * (2.0 * lz[i] + 1.0) * A2y

            dzzAy_c = (
                lz[i]
                * (lz[i] - 1.0)
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j],
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j],
                        ly[j],
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            # *** ((x-xg)(x-xk)/r^3_k + (y-yg)(y-yk)/r^3_k)dxx ***
            Aydxx_a = (
                4.0
                * exp_array[j]
                * exp_array[j]
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 3,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 2,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 2,
                        ly[j],
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 2,
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            Aydxx_b = -2.0 * exp_array[j] * (2.0 * lx[j] + 1.0) * A2y

            Aydxx_c = (
                lx[j]
                * (lx[j] - 1.0)
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] - 1,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] - 2,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] - 2,
                        ly[j],
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] - 2,
                        ly[j],
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            # *** ((x-xg)(x-xk)/r^3_k + (y-yg)(y-yk)/r^3_k)dyy ***
            Aydyy_a = (
                4.0
                * exp_array[j]
                * exp_array[j]
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 1,
                        ly[j] + 2,
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 2,
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 2,
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 2,
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            Aydyy_b = -2.0 * exp_array[j] * (2.0 * ly[j] + 1.0) * A2y

            Aydyy_c = (
                ly[j]
                * (ly[j] - 1.0)
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 1,
                        ly[j] - 2,
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] - 2,
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] - 2,
                        lz[j] + 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] - 2,
                        lz[j],
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            # *** ((x-xg)(x-xk)/r^3_k + (y-yg)(y-yk)/r^3_k)dzz ***
            Aydzz_a = (
                4.0
                * exp_array[j]
                * exp_array[j]
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j] + 2,
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 2,
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 3,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 2,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            Aydzz_b = -2.0 * exp_array[j] * (2.0 * lz[j] + 1.0) * A2y

            Aydzz_c = (
                lz[j]
                * (lz[j] - 1.0)
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j] - 2,
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] - 2,
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] - 1,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][2] - Rg[2])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] - 2,
                        0,
                        0,
                        1,
                        exp_array[i],
                        exp_array[j],
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

            intNEFy[k, i, j] = intNEFy[k, j, i] = (
                -Norm[n[i]](exp_array[i])
                * Norm[n[j]](exp_array[j])
                * 2.0
                * np.pi
                / (exp_array[i] + exp_array[j])
                * (
                    dxxAy_a
                    + dxxAy_b
                    + dxxAy_c
                    + dyyAy_a
                    + dyyAy_b
                    + dyyAy_c
                    + dzzAy_a
                    + dzzAy_b
                    + dzzAy_c
                    + Aydxx_a
                    + Aydxx_b
                    + Aydxx_c
                    + Aydyy_a
                    + Aydyy_b
                    + Aydyy_c
                    + Aydzz_a
                    + Aydzz_b
                    + Aydzz_c
                )
                * 3.0
                / 4.0
            )

            # *** dxx((x-xg)(x-xk)/r^3_k + (z-zg)(z-zk)/r^3_k) ***
            dxxAz_a = (
                4.0
                * exp_array[i]
                * exp_array[i]
                * (
                    phi.nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i] + 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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

            dxxAz_b = -2.0 * exp_array[i] * (2.0 * lx[i] + 1.0) * A2z

            dxxAz_c = (
                lx[i]
                * (lx[i] - 1.0)
                * (
                    phi.nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i] - 2,
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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

            # *** dyy((x-xg)(x-xk)/r^3_k + (z-zg)(z-zk)/r^3_k) ***
            dyyAz_a = (
                4.0
                * exp_array[i]
                * exp_array[i]
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i] + 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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

            dyyAz_b = -2.0 * exp_array[i] * (2.0 * ly[i] + 1.0) * A2z

            dyyAz_c = (
                ly[i]
                * (ly[i] - 1.0)
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i] - 2,
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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

            # *** dzz((x-xg)(x-xk)/r^3_k + (z-zg)(z-zk)/r^3_k) ***
            dzzAz_a = (
                4.0
                * exp_array[i]
                * exp_array[i]
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j],
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] + 2,
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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

            dzzAz_b = -2.0 * exp_array[i] * (2.0 * lz[i] + 1.0) * A2z

            dzzAz_c = (
                lz[i]
                * (lz[i] - 1.0)
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j] + 1,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j],
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j],
                        ly[j] + 1,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i] - 2,
                        lx[j],
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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

            # *** ((x-xg)(x-xk)/r^3_k + (z-zg)(z-zk)/r^3_k)dxx ***
            Azdxx_a = (
                4.0
                * exp_array[j]
                * exp_array[j]
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 3,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 2,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 2,
                        ly[j] + 1,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 2,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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

            Azdxx_b = -2.0 * exp_array[j] * (2.0 * lx[j] + 1.0) * A2z

            Azdxx_c = (
                lx[j]
                * (lx[j] - 1.0)
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] - 1,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] - 2,
                        ly[j],
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] - 2,
                        ly[j] + 1,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] - 2,
                        ly[j],
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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

            # *** ((x-xg)(x-xk)/r^3_k + (z-zg)(z-zk)/r^3_k)dyy ***
            Azdyy_a = (
                4.0
                * exp_array[j]
                * exp_array[j]
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 1,
                        ly[j] + 2,
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 2,
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 3,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 2,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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

            Azdyy_b = -2.0 * exp_array[j] * (2.0 * ly[j] + 1.0) * A2z

            Azdyy_c = (
                ly[j]
                * (ly[j] - 1.0)
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 1,
                        ly[j] - 2,
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] - 2,
                        lz[j],
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] - 1,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] - 2,
                        lz[j],
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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

            # *** ((x-xg)(x-xk)/r^3_k + (z-zg)(z-zk)/r^3_k)dzz ***
            Azdzz_a = (
                4.0
                * exp_array[j]
                * exp_array[j]
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j] + 2,
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 2,
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j] + 2,
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] + 2,
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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

            Azdzz_b = -2.0 * exp_array[j] * (2.0 * lz[j] + 1.0) * A2z

            Azdzz_c = (
                lz[j]
                * (lz[j] - 1.0)
                * (
                    phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + 1,
                        ly[j],
                        lz[j] - 2,
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][0] - Rg[0])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] - 2,
                        1,
                        0,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j] + 1,
                        lz[j] - 2,
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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
                    + (coord[center[j]][1] - Rg[1])
                    * phi.nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j] - 2,
                        0,
                        1,
                        0,
                        exp_array[i],
                        exp_array[j],
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

            intNEFz[k, i, j] = intNEFz[k, j, i] = (
                -Norm[n[i]](exp_array[i])
                * Norm[n[j]](exp_array[j])
                * 2.0
                * np.pi
                / (exp_array[i] + exp_array[j])
                * (
                    dxxAz_a
                    + dxxAz_b
                    + dxxAz_c
                    + dyyAz_a
                    + dyyAz_b
                    + dyyAz_c
                    + dzzAz_a
                    + dzzAz_b
                    + dzzAz_c
                    + Azdxx_a
                    + Azdxx_b
                    + Azdxx_c
                    + Azdyy_a
                    + Azdyy_b
                    + Azdyy_c
                    + Azdzz_a
                    + Azdzz_b
                    + Azdzz_c
                )
                * 3.0
                / 4.0
            )

            if output > 10:
                print("int [", i + 1, ",", j + 1, "] : ", intNEFz[k, i, j])

print("time [s] : ", time.time() - start)
