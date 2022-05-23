from numpy import exp
import numpy as np
from libh import *
from scipy.special import hyp1f1

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

# Norm = {0: phi.NS, 1: phi.NP}
intPSOx = np.zeros((2, total_nprim, total_nprim), dtype=float)
intPSOy = np.zeros((2, total_nprim, total_nprim), dtype=float)
intPSOz = np.zeros((2, total_nprim, total_nprim), dtype=float)
output = 11


for k in range(2):
    print("\n   ****Atom  ", k + 1, " ****\n")
    for i in range(total_nprim):

        for j in range(i, total_nprim):

            # pso is a combination of the NELFLD with DPVL

            # (ykdz-zkdy)/rk^3 = yk/rk^3 dz - zk/rk^3 dy =
            # V_ab^010 Dmn^1 - Vab^001 * Dkl^1  (Eq 9.931)

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

            zdy = 2.0 * exp[j] * nuclear_attraction(
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
            ) - ly[j] * nuclear_attraction(
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

            intPSOx[k, i, j] = intPSOx[k, j, i] = (
                -Norm[lx[i] + ly[i] + lz[i]](exp[i])
                * Norm[lx[j] + ly[j] + lz[j]](exp[j])
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
                * (ydz - zdy)
            )

            # (zkdx-xkdz)/rk^3 = zk/rk^3 dx - xk/rk^3 dz =
            # V_ab^001 Dij^1 - Vab^100 * Dmn^1  (Eq 9.931)

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

            intPSOy[k, i, j] = intPSOy[k, j, i] = (
                -Norm[lx[i] + ly[i] + lz[i]](exp[i])
                * Norm[lx[j] + ly[j] + lz[j]](exp[j])
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
                * (zdx - xdz)
            )

            # (xkdy-ykdz)/rk^3 = xk/rk^3 dy - yk/rk^3 dx =
            # V_ab^100 Dkl^1 - Vab^010 * Dij^1  (Eq 9.931)

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

            intPSOz[k, i, j] = intPSOz[k, j, i] = (
                -Norm[lx[i] + ly[i] + lz[i]](exp[i])
                * Norm[lx[j] + ly[j] + lz[j]](exp[j])
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
                * (xdy - ydx)
            )

            if output > 10 and np.abs(intPSOy[k, i, j]) > 1e-2:
                print("int [", i + 1, ",", j + 1, "] : ", intPSOy[k, i, j])

print("time [s] : ", time() - start)
