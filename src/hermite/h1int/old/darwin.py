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
intDwx = np.zeros((total_nprim, total_nprim), dtype=float)
output = 11

c = 137.0359998
a2 = 1.0 / (c ** 2)
darwin_const = np.pi * a2 / 2.0


def gaussian_mult(
    i, k, m, j, l, n, Ax, Ay, Az, Bx, By, Bz, alpha, beta, Kx, Ky, Kz
):
    # <phi|delta(rk)|phi> = <phi|delta(x-Kx,y-Ky,z-Kz)|phi> =
    # (xk-Ax)^i(yk-Ay)^k(zk-Az)^mexp(-alpha*(rk-A)^2) *
    # (xk-Bx)^j(yk-By)^l(zk-Bz)^nexp(-alpha*(rk-B)^2)

    gij = exp(-alpha * (Kx - Ax) ** 2) * exp(-beta * (Kx - Bx) ** 2)
    #    gij = np.power((Kx - Ax), i) * np.power((Kx - Bx), j) * gij
    # !TODO review this change in order of (Ay, By, Ky) y (Ax, Bx, Kx)
    gij = np.power((Ax - Kx), i) * np.power((Bx - Kx), j) * gij

    gkl = exp(-alpha * (Ky - Ay) ** 2) * exp(-beta * (Ky - By) ** 2)
    #    gkl = np.power((Ky - Ay), k) * np.power((Ky - By), l) * gkl
    gkl = np.power((Ay - Ky), k) * np.power((By - Ky), l) * gkl

    gmn = exp(-alpha * (Kz - Az) ** 2) * exp(-beta * (Kz - Bz) ** 2)
    #    gmn = np.power((Az - Kz), m) * np.power((Bz - Kz), n) * gmn
    # ! To reproduce sign of DALTON like with Fc, it is change x - Ax by Ax - x
    gmn = np.power((Az - Kz), m) * np.power((Bz - Kz), n) * gmn

    return gij * gkl * gmn


for k in range(2):
    # print("\n   ****Atom  ", k + 1, " ****\n")
    for i in range(total_nprim):

        for j in range(i, total_nprim):

            dw = gaussian_mult(
                lx[i],
                ly[i],
                lz[i],
                lx[j],
                ly[j],
                lz[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                exp_array[i],
                exp_array[j],
                coord[k][0],
                coord[k][1],
                coord[k][2],
            )

            intDwx[i, j] += (
                darwin_const
                * Norm[n[i]](exp_array[i])
                * Norm[n[j]](exp_array[j])
                * dw
            )
            intDwx[j, i] = intDwx[i, j]

            if output > 10 and k == 1:
                print("int [", j + 1, ",", i + 1, "] : ", intDwx[j, i])

print(" time [s]: ", -start + time.time())
