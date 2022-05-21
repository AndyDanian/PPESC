from numpy import exp
import numpy as np
#import time
from libh import *

start = time()

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

lx = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]
ly = [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
lz = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1]

center = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1]
total_nprim = 18

coord = [[0.0, 0.0, 0.0586476414], [0.0, 0.0, 1.4045523587]]

intLx = np.zeros((total_nprim, total_nprim), dtype=float)
intLy = np.zeros((total_nprim, total_nprim), dtype=float)
intLz = np.zeros((total_nprim, total_nprim), dtype=float)
output = 11

# !Note: Rgaugeo is put on in the origen
Rg = [0.0, 0.0, 1.404552358700]
for i in range(total_nprim):

    for j in range(i, total_nprim):

        Px = (
            exp_array[i] * coord[center[i]][0]
            + exp_array[j] * coord[center[j]][0]
        )
        Px = Px / (exp_array[i] + exp_array[j])

        Py = (
            exp_array[i] * coord[center[i]][1]
            + exp_array[j] * coord[center[j]][1]
        )
        Py = Py / (exp_array[i] + exp_array[j])

        Pz = (
            exp_array[i] * coord[center[i]][2]
            + exp_array[j] * coord[center[j]][2]
        )
        Pz = Pz / (exp_array[i] + exp_array[j])

        Xpg = Px - Rg[0]
        Ypg = Py - Rg[1]
        Zpg = Pz - Rg[2]
        print("ygaugeo zgaugeo ", Ypg, Zpg)


        sij = E(
            lx[i],
            lx[j],
            0,
            coord[center[i]][0] - coord[center[j]][0],
            exp_array[i],
            exp_array[j],
        )

        skl = E(
            ly[i],
            ly[j],
            0,
            coord[center[i]][1] - coord[center[j]][1],
            exp_array[i],
            exp_array[j],
        )

        smn = E(
            lz[i],
            lz[j],
            0,
            coord[center[i]][2] - coord[center[j]][2],
            exp_array[i],
            exp_array[j],
        )

        # Real{Lx} = <phi|(ygdz - zgdy)|phi> =
        #  [(ygSkl^0)(dzSmn^0) - (zgSmn^0)(dySkl^0)]Sij^0 =
        #  [Skl^1Dmn^1 - Smn^1Dkl^1]Sij^0 =
        #  [(E_1^kl + YpgE0^kl)(2bE_0^mn+1-l2bE0^mn-1) -
        #   (E_1^mn + ZpgE0^mn)(2bE_0^kl+1-l2bE0^kl-1)]Sij^0
        ygaugeo = E(
            ly[i],
            ly[j],
            1,
            coord[center[i]][1] - coord[center[j]][1],
            exp_array[i],
            exp_array[j],
        )

        py = 2.0 * exp_array[j] * E(
            ly[i],
            ly[j] + 1,
            0,
            coord[center[i]][1] - coord[center[j]][1],
            exp_array[i],
            exp_array[j],
        ) - ly[j] * E(
            ly[i],
            ly[j] - 1,
            0,
            coord[center[i]][1] - coord[center[j]][1],
            exp_array[i],
            exp_array[j],
        )

        zgaugeo = E(
            lz[i],
            lz[j],
            1,
            coord[center[i]][2] - coord[center[j]][2],
            exp_array[i],
            exp_array[j],
        )

        pz = 2.0 * exp_array[j] * E(
            lz[i],
            lz[j] + 1,
            0,
            coord[center[i]][2] - coord[center[j]][2],
            exp_array[i],
            exp_array[j],
        ) - lz[j] * E(
            lz[i],
            lz[j] - 1,
            0,
            coord[center[i]][2] - coord[center[j]][2],
            exp_array[i],
            exp_array[j],
        )

        intLx[i, j] = (
            Norm[lx[i] + ly[i] + lz[i]](exp_array[i])
            * Norm[lx[j] + ly[j] + lz[j]](exp_array[j])
            * ((ygaugeo + Ypg * skl) * pz - (zgaugeo + Zpg * smn) * py)
            * sij
            * np.power(np.pi / (exp_array[i] + exp_array[j]), 1.5)
        )
        intLx[j, i] = -1.0 * intLx[i, j]

        # Real{Ly} = <phi|(zgdx - xgdz)|phi> =
        #  [(zgSmn^0)(dxSij^0) - (xgSij^0)(dzSmn^0)]Skl^0 =
        #  [Smn^1Dij^1 - Sij^1Dmn^1]Skl^0 =
        #  [(E_1^mn + ZpgE0^mn)(2bE_0^ij+1-l2bE0^ij-1) -
        #   (E_1^ij + XpgE0^ij)(2bE_0^mn+1-l2bE0^mn-1)]Skl^0

        xgaugeo = E(
            lx[i],
            lx[j],
            1,
            coord[center[i]][0] - coord[center[j]][0],
            exp_array[i],
            exp_array[j],
        )

        px = 2.0 * exp_array[j] * E(
            lx[i],
            lx[j] + 1,
            0,
            coord[center[i]][0] - coord[center[j]][0],
            exp_array[i],
            exp_array[j],
        ) - lx[j] * E(
            lx[i],
            lx[j] - 1,
            0,
            coord[center[i]][0] - coord[center[j]][0],
            exp_array[i],
            exp_array[j],
        )

        intLy[i, j] = (
            Norm[lx[i] + ly[i] + lz[i]](exp_array[i])
            * Norm[lx[j] + ly[j] + lz[j]](exp_array[j])
            * ((zgaugeo + Zpg * smn) * px - (xgaugeo + Xpg * sij) * pz)
            * skl
            * np.power(np.pi / (exp_array[i] + exp_array[j]), 1.5)
        )
        intLy[j, i] = -1.0 * intLy[i, j]

        # Real{Lz} = <phi|(xgdy - ygdx)|phi> =
        #  [(xgSij^0)(dySkl^0) - (ygSkl^0)(dxSij^0)]Smn^0 =
        #  [Sij^1Dkl^1 - Skl^1Dij^1]Smn^0 =
        #  [(E_1^ij + XpgE0^ij)(2bE_0^kl+1-l2bE0^kl-1) -
        #   (E_1^kl + YpgE0^kl)(2bE_0^ij+1-l2bE0^ij-1)]Smn^0

        intLz[i, j] = (
            Norm[lx[i] + ly[i] + lz[i]](exp_array[i])
            * Norm[lx[j] + ly[j] + lz[j]](exp_array[j])
            * ((xgaugeo + Xpg * sij) * py - (ygaugeo + Ypg * skl) * px)
            * smn
            * np.power(np.pi / (exp_array[i] + exp_array[j]), 1.5)
        )
        intLz[j, i] = -1.0 * intLz[i, j]
        # ! Nota: Revisar que tipo de matriz es, summetric, antisymmetric or Cuadrada
        if output > 10 and np.abs(intLx[i, j]) > 1e-2:
            print("int [", i + 1, ",", j + 1, "] : ", intLx[i, j])

print(" time [s]: ", -start + time())
