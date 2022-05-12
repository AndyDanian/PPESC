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
intDx = np.zeros((total_nprim, total_nprim), dtype=float)
intDy = np.zeros((total_nprim, total_nprim), dtype=float)
intDz = np.zeros((total_nprim, total_nprim), dtype=float)
output = 11

# !Note: rc is put on in the origen
Rdipole = [0.0, 0.0, 0.0]

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
        # Eq 9.5.43 Helgaker
        # <phi|x_c|phi> = <phi|x_p + Xpc|phi> =
        # Sij^1Skl^0Smn^0=(E1^ij + Xpc*E0^ij)E0^klE0^mn(pi/p)^1.5

        Px = (
            exp_array[i] * coord[center[i]][0]
            + exp_array[j] * coord[center[j]][0]
        )
        Px = Px / (exp_array[i] + exp_array[j])
        Xpk = Px - Rdipole[0]

        xdipole = phi.E(
            lx[i],
            lx[j],
            1,
            coord[center[i]][0] - coord[center[j]][0],
            exp_array[i],
            exp_array[j],
        )

        intDx[i, j] = intDx[j, i] = (
            Norm[n[i]](exp_array[i])
            * Norm[n[j]](exp_array[j])
            * (xdipole + Xpk * sij)
            * skl
            * smn
            * np.power(np.pi / (exp_array[i] + exp_array[j]), 1.5)
        )

        # <phi|y_c|phi> = <phi|y_p + Ypc|phi> =
        # Sij^0Skl^1Smn^0=E0^ij(E1^kl + XpcE0^kl)E0^mn(pi/p)^1.5
        Py = (
            exp_array[i] * coord[center[i]][1]
            + exp_array[j] * coord[center[j]][1]
        )
        Py = Py / (exp_array[i] + exp_array[j])
        Ypk = Py - Rdipole[1]

        ydipole = phi.E(
            ly[i],
            ly[j],
            1,
            coord[center[i]][1] - coord[center[j]][1],
            exp_array[i],
            exp_array[j],
        )

        intDy[i, j] = intDy[j, i] = (
            Norm[n[i]](exp_array[i])
            * Norm[n[j]](exp_array[j])
            * (ydipole + Ypk * skl)
            * sij
            * smn
            * np.power(np.pi / (exp_array[i] + exp_array[j]), 1.5)
        )

        # <phi|z_c|phi> = <phi|z_p + Zpc|phi> =
        # Sij^0Skl^0Smn^1=E0^ijE0^kl(E1^mn + XpcE0^mn)(pi/p)^1.5
        Pz = (
            exp_array[i] * coord[center[i]][2]
            + exp_array[j] * coord[center[j]][2]
        )
        Pz = Pz / (exp_array[i] + exp_array[j])
        Zpk = Pz - Rdipole[2]

        dz = phi.E(
            lz[i],
            lz[j],
            1,
            coord[center[i]][2] - coord[center[j]][2],
            exp_array[i],
            exp_array[j],
        )

        intDz[i, j] = intDz[j, i] = (
            Norm[n[i]](exp_array[i])
            * Norm[n[j]](exp_array[j])
            * (dz + Zpk * smn)
            * sij
            * skl
            * np.power(np.pi / (exp_array[i] + exp_array[j]), 1.5)
        )

        if output > 10 and np.abs(intDx[i, j]) > 1e-2:
            print("int [", i + 1, ",", j + 1, "] : ", intDx[i, j])

print(" time [s]: ", -start + time.time())
