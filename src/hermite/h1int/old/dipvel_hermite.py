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
intPx = np.zeros((total_nprim, total_nprim), dtype=float)
intPy = np.zeros((total_nprim, total_nprim), dtype=float)
intPz = np.zeros((total_nprim, total_nprim), dtype=float)
output = 11

for i in range(total_nprim):

    for j in range(i, total_nprim):

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

        # Horizontal Reccurence

        # int phi d/dx phi dt = Dij^1Skl^0Smn^0 = (2*b*Sij+1^0 - j*Sij-1^0)Skl^0Smn^0
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
        intPx[i, j] = (
            Norm[n[i]](exp_array[i])
            * Norm[n[j]](exp_array[j])
            * px
            * skl
            * smn
            * np.power(np.pi / (exp_array[i] + exp_array[j]), 1.5)
        )
        intPx[j, i] = -1.0 * intPx[i, j]

        # int phi d/dy phi dt = Sij^1Dkl^1Smn^0 = Sij^0(2*b*Skl+1^0 - j*Skl-1^0)Smn^0
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
        intPy[i, j] = (
            Norm[n[i]](exp_array[i])
            * Norm[n[j]](exp_array[j])
            * py
            * sij
            * smn
            * np.power(np.pi / (exp_array[i] + exp_array[j]), 1.5)
        )
        intPy[j, i] = -1.0 * intPy[i, j]

        # int phi d/dz phi dt = Sij^0Skl^0Dmn^l = Sij^0Skl^0(2*b*Smn+1^1 - j*Smn-1^1)
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
        intPz[i, j] = (
            Norm[n[i]](exp_array[i])
            * Norm[n[j]](exp_array[j])
            * pz
            * sij
            * skl
            * np.power(np.pi / (exp_array[i] + exp_array[j]), 1.5)
        )
        intPz[j, i] = -1.0 * intPz[i, j]

        # Pi Matrices are symmetric, but one half is the negative of the another
        if output > 10 and np.abs(intPx[i, j]) > 1e-3:
            print("int [", i + 1, ",", j + 1, "] : ", intPx[i, j])

print(" time [s]: ", -start + time.time())
