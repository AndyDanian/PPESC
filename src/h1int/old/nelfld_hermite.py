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


for k in range(2):
    print("\n   ****Atom  ", k + 1, " ****\n")
    for i in range(total_nprim):

        for j in range(i, total_nprim):

            # Nuclear Electric Field Gradient
            # <phi|xk/rk^3||phi> --> Eqs 9.932, 9.9.18-20

            # *** x
            nef = phi.nuclear_attraction(
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

            intNEFx[k, i, j] = intNEFx[k, j, i] = (
                -Norm[n[i]](exp_array[i])
                * Norm[n[j]](exp_array[j])
                * 2.0
                * np.pi
                / (exp_array[i] + exp_array[j])
                * nef
            )

            # *** y
            nef = phi.nuclear_attraction(
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

            intNEFy[k, i, j] = intNEFy[k, j, i] = (
                -Norm[n[i]](exp_array[i])
                * Norm[n[j]](exp_array[j])
                * 2.0
                * np.pi
                / (exp_array[i] + exp_array[j])
                * nef
            )

            # *** z
            nef = phi.nuclear_attraction(
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

            intNEFz[k, i, j] = intNEFz[k, j, i] = (
                -Norm[n[i]](exp_array[i])
                * Norm[n[j]](exp_array[j])
                * 2.0
                * np.pi
                / (exp_array[i] + exp_array[j])
                * nef
            )

            if output > 10:
                print("int [", i + 1, ",", j + 1, "] : ", intNEFz[k, i, j])

print("time [s] : ", time.time() - start)
