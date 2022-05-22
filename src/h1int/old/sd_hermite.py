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
intSDx = np.zeros((2, total_nprim, total_nprim), dtype=float)
intSDy = np.zeros((2, total_nprim, total_nprim), dtype=float)
intSDz = np.zeros((2, total_nprim, total_nprim), dtype=float)
output = 11

# ! Buscar ref gfactor
gfactor = 2.0023193134
const_sd = gfactor / 2.0 * 1 / 3.0


for k in range(2):
    print("\n   ****Atom  ", k + 1, " ****\n")
    for i in range(total_nprim):

        for j in range(i, total_nprim):

            # SD = <phi|(3r_krk^T-r_k^2)/rk^5|phi>
            # !Note: is add the constant 1/3. to reproduce of DALTON value

            # *** x
            x2r5 = phi.nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j],
                ly[j],
                lz[j],
                2,
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
            y2r5 = phi.nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j],
                ly[j],
                lz[j],
                0,
                2,
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
            z2r5 = phi.nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j],
                ly[j],
                lz[j],
                0,
                0,
                2,
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

            intSDx[k, i, j] = intSDx[k, j, i] = (
                const_sd
                * Norm[n[i]](exp_array[i])
                * Norm[n[j]](exp_array[j])
                * 2.0
                * np.pi
                / (exp_array[i] + exp_array[j])
                * (2.0 * x2r5 - y2r5 - z2r5)
            )

            # *** y
            intSDy[k, i, j] = intSDy[k, j, i] = (
                const_sd
                * Norm[n[i]](exp_array[i])
                * Norm[n[j]](exp_array[j])
                * 2.0
                * np.pi
                / (exp_array[i] + exp_array[j])
                * (2.0 * y2r5 - x2r5 - z2r5)
            )

            # *** z
            intSDz[k, i, j] = intSDz[k, j, i] = (
                const_sd
                * Norm[n[i]](exp_array[i])
                * Norm[n[j]](exp_array[j])
                * 2.0
                * np.pi
                / (exp_array[i] + exp_array[j])
                * (2.0 * z2r5 - x2r5 - y2r5)
            )

            if output > 10 and np.abs(intSDz[k, i, j]) > 1e-2:
                print("int [", i + 1, ",", j + 1, "] : ", intSDz[k, i, j])

print("time [s] : ", time.time() - start)
