from numpy import exp
import numpy as np
import time
import phi
from scipy.special import hyp1f1

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
intPSOx = np.zeros((2, total_nprim, total_nprim), dtype=float)
intPSOy = np.zeros((2, total_nprim, total_nprim), dtype=float)
intPSOz = np.zeros((2, total_nprim, total_nprim), dtype=float)
output = 11

Rg = [0.0, 0.0, 1.404552358700]


for k in range(2):
    print("\n   ****Atom  ", k + 1, " ****\n")
    for i in range(total_nprim):

        for j in range(total_nprim):

            # (ygdz-zgdy)(ykdz-zkdy)/rk^3 =
            # (ygdz) yk/rk^3 dz - (ygdz) zk/rk^3 dy - (zgdy) yk/rk^3 dz + (zgdy) zk/rk^3 dy=
            # (yk + ky - yg) Dmn^1 V_ab^010 Dmn^1 - (yk + ky - yg) Dmn^1 Vab^001 Dkl^1
            # (zk + kz - zg) Dkl^1 V_ab^010 Dmn^1 - (zk + kz - zg) Dkl^1 Vab^001 Dkl^1 (Eq 9.931)

            # * yAdz yk/r^3k dz
            l0xlkx_a = (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i] + 1,
                    lx[j],
                    ly[j],
                    lz[j] + 1,
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
                - 2.0
                * exp_array[i]
                * lz[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i] + 1,
                    lx[j],
                    ly[j],
                    lz[j] - 1,
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
                - 2.0
                * exp_array[j]
                * lz[i]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i] - 1,
                    lx[j],
                    ly[j],
                    lz[j] + 1,
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
                + lz[i]
                * lz[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i] - 1,
                    lx[j],
                    ly[j],
                    lz[j] - 1,
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

            # * (Ay-Gy)dz yk/r^3k dz
            l0xlkx_b = (coord[center[i]][1] - Rg[1]) * (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i] + 1,
                    lx[j],
                    ly[j],
                    lz[j] + 1,
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
                - 2.0
                * exp_array[i]
                * lz[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i] + 1,
                    lx[j],
                    ly[j],
                    lz[j] - 1,
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
                - 2.0
                * exp_array[j]
                * lz[i]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i] - 1,
                    lx[j],
                    ly[j],
                    lz[j] + 1,
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
                + lz[i]
                * lz[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i] - 1,
                    lx[j],
                    ly[j],
                    lz[j] - 1,
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

            # * zAdy yk/r^3k dz
            l0xlkx_c = (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i] + 1,
                    lx[j],
                    ly[j],
                    lz[j] + 1,
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
                - 2.0
                * exp_array[i]
                * lz[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i] + 1,
                    lx[j],
                    ly[j],
                    lz[j] - 1,
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
                - 2.0
                * exp_array[j]
                * ly[i]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] - 1,
                    lz[i] + 1,
                    lx[j],
                    ly[j],
                    lz[j] + 1,
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
                + ly[i]
                * lz[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] - 1,
                    lz[i] + 1,
                    lx[j],
                    ly[j],
                    lz[j] - 1,
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

            # * (Az-Gz)dy yk/r^3k dz
            l0xlkx_d = (coord[center[i]][2] - Rg[2]) * (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i],
                    lx[j],
                    ly[j],
                    lz[j] + 1,
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
                - 2.0
                * exp_array[i]
                * lz[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i],
                    lx[j],
                    ly[j],
                    lz[j] - 1,
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
                - 2.0
                * exp_array[j]
                * ly[i]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] - 1,
                    lz[i],
                    lx[j],
                    ly[j],
                    lz[j] + 1,
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
                + ly[i]
                * lz[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] - 1,
                    lz[i],
                    lx[j],
                    ly[j],
                    lz[j] - 1,
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

            # * yAdz zk/r^3k dy
            l0xlkx_e = (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i] + 1,
                    lx[j],
                    ly[j] + 1,
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
                - 2.0
                * exp_array[i]
                * ly[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i] + 1,
                    lx[j],
                    ly[j] - 1,
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
                - 2.0
                * exp_array[j]
                * lz[i]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i] - 1,
                    lx[j],
                    ly[j] + 1,
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
                + lz[i]
                * ly[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i] - 1,
                    lx[j],
                    ly[j] - 1,
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

            # * (Ay-Gy)dz zk/r^3k dy
            l0xlkx_f = (coord[center[i]][1] - Rg[1]) * (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i] + 1,
                    lx[j],
                    ly[j] + 1,
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
                - 2.0
                * exp_array[i]
                * ly[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i] + 1,
                    lx[j],
                    ly[j] - 1,
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
                - 2.0
                * exp_array[j]
                * lz[i]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i] - 1,
                    lx[j],
                    ly[j] + 1,
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
                + lz[i]
                * ly[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i] - 1,
                    lx[j],
                    ly[j] - 1,
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

            # * zAdy zk/r^3k dy
            l0xlkx_g = (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i] + 1,
                    lx[j],
                    ly[j] + 1,
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
                - 2.0
                * exp_array[i]
                * ly[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i] + 1,
                    lx[j],
                    ly[j] - 1,
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
                - 2.0
                * exp_array[j]
                * ly[i]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] - 1,
                    lz[i] + 1,
                    lx[j],
                    ly[j] + 1,
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
                + ly[i]
                * ly[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] - 1,
                    lz[i] + 1,
                    lx[j],
                    ly[j] - 1,
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

            # * (Ay-Gy)dz zk/r^3k dy
            l0xlkx_h = (coord[center[i]][2] - Rg[2]) * (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i],
                    lx[j],
                    ly[j] + 1,
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
                - 2.0
                * exp_array[i]
                * ly[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i],
                    lx[j],
                    ly[j] - 1,
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
                - 2.0
                * exp_array[j]
                * ly[i]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] - 1,
                    lz[i],
                    lx[j],
                    ly[j] + 1,
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
                + ly[i]
                * ly[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] - 1,
                    lz[i],
                    lx[j],
                    ly[j] - 1,
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

            intPSOx[k, j, i] = (
                -Norm[n[i]](exp_array[i])
                * Norm[n[j]](exp_array[j])
                * 2.0
                * np.pi
                / (exp_array[i] + exp_array[j])
                * (
                    l0xlkx_a
                    + l0xlkx_b
                    - l0xlkx_c
                    - l0xlkx_d
                    - l0xlkx_e
                    - l0xlkx_f
                    + l0xlkx_g
                    + l0xlkx_h
                )
            )

            # (zgdx-xgdz)(zkdx-xkdz)/rk^3 =
            # (zgdx) zk/rk^3 dx - (zgdx) xk/rk^3 dz - (xgdz) zk/rk^3 dx + (xgdz) xk/rk^3 dz=
            # (zk + kz - zg) Dij^1 V_ab^001 Dij^1 - (zk + kz - zg) Dij^1 Vab^100 Dmn^1
            # (xk + xz - xg) Dmn^1 V_ab^001 Dij^1 - (xk + kx - xg) Dmn^1 Vab^100 Dmn^1

            # * zA dx zk/r^3k dx
            l0ylky_a = (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i] + 1,
                    lx[j] + 1,
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
                - 2.0
                * exp_array[i]
                * lx[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i] + 1,
                    lx[j] - 1,
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
                - 2.0
                * exp_array[j]
                * lx[i]
                * phi.nuclear_attraction(
                    lx[i] - 1,
                    ly[i],
                    lz[i] + 1,
                    lx[j] + 1,
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
                + lx[i]
                * lx[j]
                * phi.nuclear_attraction(
                    lx[i] - 1,
                    ly[i],
                    lz[i] + 1,
                    lx[j] - 1,
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

            # * (Az-Gz)dx zk/r^3k dx
            l0ylky_b = (coord[center[i]][2] - Rg[2]) * (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i],
                    lx[j] + 1,
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
                - 2.0
                * exp_array[i]
                * lx[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i],
                    lx[j] - 1,
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
                - 2.0
                * exp_array[j]
                * lx[i]
                * phi.nuclear_attraction(
                    lx[i] - 1,
                    ly[i],
                    lz[i],
                    lx[j] + 1,
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
                + lx[i]
                * lx[j]
                * phi.nuclear_attraction(
                    lx[i] - 1,
                    ly[i],
                    lz[i],
                    lx[j] - 1,
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

            # * zA dx xk/r^3k dz
            l0ylky_c = (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i] + 1,
                    lx[j],
                    ly[j],
                    lz[j] + 1,
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
                - 2.0
                * exp_array[i]
                * lz[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i] + 1,
                    lx[j],
                    ly[j],
                    lz[j] - 1,
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
                - 2.0
                * exp_array[j]
                * lx[i]
                * phi.nuclear_attraction(
                    lx[i] - 1,
                    ly[i],
                    lz[i] + 1,
                    lx[j],
                    ly[j],
                    lz[j] + 1,
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
                + lx[i]
                * lz[j]
                * phi.nuclear_attraction(
                    lx[i] - 1,
                    ly[i],
                    lz[i] + 1,
                    lx[j],
                    ly[j],
                    lz[j] - 1,
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
            )

            # * (Az-Gz) dx xk/r^3k dz
            l0ylky_d = (coord[center[i]][2] - Rg[2]) * (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j],
                    lz[j] + 1,
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
                - 2.0
                * exp_array[i]
                * lz[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j],
                    lz[j] - 1,
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
                - 2.0
                * exp_array[j]
                * lx[i]
                * phi.nuclear_attraction(
                    lx[i] - 1,
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j],
                    lz[j] + 1,
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
                + lx[i]
                * lz[j]
                * phi.nuclear_attraction(
                    lx[i] - 1,
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j],
                    lz[j] - 1,
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
            )

            # * xA dz zk/r^3k dx
            l0ylky_e = (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i] + 1,
                    lx[j] + 1,
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
                - 2.0
                * exp_array[i]
                * lx[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i] + 1,
                    lx[j] - 1,
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
                - 2.0
                * exp_array[j]
                * lz[i]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i] - 1,
                    lx[j] + 1,
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
                + lz[i]
                * lx[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i] - 1,
                    lx[j] - 1,
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

            # * (Ax-Gx) dz zk/r^3k dx
            l0ylky_f = (coord[center[i]][0] - Rg[0]) * (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i] + 1,
                    lx[j] + 1,
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
                - 2.0
                * exp_array[i]
                * lx[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i] + 1,
                    lx[j] - 1,
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
                - 2.0
                * exp_array[j]
                * lz[i]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i] - 1,
                    lx[j] + 1,
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
                + lz[i]
                * lx[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i] - 1,
                    lx[j] - 1,
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

            # * xA dz xk/r^3k dz
            l0ylky_g = (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i] + 1,
                    lx[j],
                    ly[j],
                    lz[j] + 1,
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
                - 2.0
                * exp_array[i]
                * lz[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i] + 1,
                    lx[j],
                    ly[j],
                    lz[j] - 1,
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
                - 2.0
                * exp_array[j]
                * lz[i]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i] - 1,
                    lx[j],
                    ly[j],
                    lz[j] + 1,
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
                + lz[i]
                * lz[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i] - 1,
                    lx[j],
                    ly[j],
                    lz[j] - 1,
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
            )

            # * (Ax-Gx) dz xk/r^3k dz
            l0ylky_h = (coord[center[i]][0] - Rg[0]) * (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i] + 1,
                    lx[j],
                    ly[j],
                    lz[j] + 1,
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
                - 2.0
                * exp_array[i]
                * lz[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i] + 1,
                    lx[j],
                    ly[j],
                    lz[j] - 1,
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
                - 2.0
                * exp_array[j]
                * lz[i]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i] - 1,
                    lx[j],
                    ly[j],
                    lz[j] + 1,
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
                + lz[i]
                * lz[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i] - 1,
                    lx[j],
                    ly[j],
                    lz[j] - 1,
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
            )

            intPSOy[k, j, i] = (
                -Norm[n[i]](exp_array[i])
                * Norm[n[j]](exp_array[j])
                * 2.0
                * np.pi
                / (exp_array[i] + exp_array[j])
                * (
                    l0ylky_a
                    + l0ylky_b
                    - l0ylky_c
                    - l0ylky_d
                    - l0ylky_e
                    - l0ylky_f
                    + l0ylky_g
                    + l0ylky_h
                )
            )

            # (xgdy-ygdx)(xkdy-ykdx)/rk^3 =
            # (xgdy) xk/rk^3 dy - (xgdy) yk/rk^3 dx - (ygdx) xk/rk^3 dy + (ygdx) yk/rk^3 dx=
            # (xk + kx - xg) Dkl^1 V_ab^100 Dkl^1 - (xk + kx - xg) Dkl^1 Vab^010 Dij^1
            # (yk + ky - yg) Dij^1 V_ab^100 Dkl^1 - (yk + ky - yg) Dij^1 Vab^010 Dij^1

            # * xA dy xk/r^3k dy
            l0zlkz_a = (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i] + 1,
                    lz[i],
                    lx[j],
                    ly[j] + 1,
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
                - 2.0
                * exp_array[i]
                * ly[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i] + 1,
                    lz[i],
                    lx[j],
                    ly[j] - 1,
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
                - 2.0
                * exp_array[j]
                * ly[i]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i] - 1,
                    lz[i],
                    lx[j],
                    ly[j] + 1,
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
                + ly[i]
                * ly[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i] - 1,
                    lz[i],
                    lx[j],
                    ly[j] - 1,
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
            )

            # * (Az-Gz)dx zk/r^3k dx
            l0zlkz_b = (coord[center[i]][0] - Rg[0]) * (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i],
                    lx[j],
                    ly[j] + 1,
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
                - 2.0
                * exp_array[i]
                * ly[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i],
                    lx[j],
                    ly[j] - 1,
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
                - 2.0
                * exp_array[j]
                * ly[i]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] - 1,
                    lz[i],
                    lx[j],
                    ly[j] - 1,
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
                + ly[i]
                * ly[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] - 1,
                    lz[i],
                    lx[j],
                    ly[j] - 1,
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
            )

            # * xA dy yk/r^3k dx
            l0zlkz_c = (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i] + 1,
                    lz[i],
                    lx[j] + 1,
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
                - 2.0
                * exp_array[i]
                * lx[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i] + 1,
                    lz[i],
                    lx[j] - 1,
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
                - 2.0
                * exp_array[j]
                * ly[i]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i] - 1,
                    lz[i],
                    lx[j] + 1,
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
                + ly[i]
                * lx[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i] - 1,
                    lz[i],
                    lx[j] - 1,
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

            # * (Ax-Gx) dy yk/r^3k dx
            l0zlkz_d = (coord[center[i]][0] - Rg[0]) * (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i],
                    lx[j] + 1,
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
                - 2.0
                * exp_array[i]
                * lx[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] + 1,
                    lz[i],
                    lx[j] - 1,
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
                - 2.0
                * exp_array[j]
                * ly[i]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] - 1,
                    lz[i],
                    lx[j] + 1,
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
                + ly[i]
                * lx[j]
                * phi.nuclear_attraction(
                    lx[i],
                    ly[i] - 1,
                    lz[i],
                    lx[j] - 1,
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

            # * yA dx xk/r^3k dy
            l0zlkz_e = (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i] + 1,
                    lz[i],
                    lx[j],
                    ly[j] + 1,
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
                - 2.0
                * exp_array[i]
                * ly[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i] + 1,
                    lz[i],
                    lx[j],
                    ly[j] - 1,
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
                - 2.0
                * exp_array[j]
                * lx[i]
                * phi.nuclear_attraction(
                    lx[i] - 1,
                    ly[i] + 1,
                    lz[i],
                    lx[j],
                    ly[j] + 1,
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
                + lx[i]
                * ly[j]
                * phi.nuclear_attraction(
                    lx[i] - 1,
                    ly[i] + 1,
                    lz[i],
                    lx[j],
                    ly[j] - 1,
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
            )

            # * (Ay-Gy) dx xk/r^3k dy
            l0zlkz_f = (coord[center[i]][1] - Rg[1]) * (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j] + 1,
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
                - 2.0
                * exp_array[i]
                * ly[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j] - 1,
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
                - 2.0
                * exp_array[j]
                * lx[i]
                * phi.nuclear_attraction(
                    lx[i] - 1,
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j] + 1,
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
                + lx[i]
                * ly[j]
                * phi.nuclear_attraction(
                    lx[i] - 1,
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j] - 1,
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
            )

            # * yA dx yk/r^3k dx
            l0zlkz_g = (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i] + 1,
                    lz[i],
                    lx[j] + 1,
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
                - 2.0
                * exp_array[i]
                * lx[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i] + 1,
                    lz[i],
                    lx[j] - 1,
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
                - 2.0
                * exp_array[j]
                * lx[i]
                * phi.nuclear_attraction(
                    lx[i] - 1,
                    ly[i] + 1,
                    lz[i],
                    lx[j] + 1,
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
                + lx[i]
                * lx[j]
                * phi.nuclear_attraction(
                    lx[i] - 1,
                    ly[i] + 1,
                    lz[i],
                    lx[j] - 1,
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

            # * (Ay-Gy) dx yk/r^3k dx
            l0zlkz_h = (coord[center[i]][1] - Rg[1]) * (
                4.0
                * exp_array[i]
                * exp_array[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i],
                    lx[j] + 1,
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
                - 2.0
                * exp_array[i]
                * lx[j]
                * phi.nuclear_attraction(
                    lx[i] + 1,
                    ly[i],
                    lz[i],
                    lx[j] - 1,
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
                - 2.0
                * exp_array[j]
                * lx[i]
                * phi.nuclear_attraction(
                    lx[i] - 1,
                    ly[i],
                    lz[i],
                    lx[j] + 1,
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
                + lx[i]
                * lx[j]
                * phi.nuclear_attraction(
                    lx[i] - 1,
                    ly[i],
                    lz[i],
                    lx[j] - 1,
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

            intPSOz[k, j, i] = (
                -Norm[n[i]](exp_array[i])
                * Norm[n[j]](exp_array[j])
                * 2.0
                * np.pi
                / (exp_array[i] + exp_array[j])
                * (
                    l0zlkz_a
                    + l0zlkz_b
                    - l0zlkz_c
                    - l0zlkz_d
                    - l0zlkz_e
                    - l0zlkz_f
                    + l0zlkz_g
                    + l0zlkz_h
                )
            )

            if output > 10 and np.abs(intPSOz[k, j, i]) > 1e-2:
                print("int [", j + 1, ",", i + 1, "] : ", intPSOz[k, j, i])

print("time [s] : ", time.time() - start)
