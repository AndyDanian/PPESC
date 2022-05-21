from libh import *

# from numpy import exp
# import numpy as np
# import time
# import phi

# start = time.time()

# # 6-311++G**
# exp = [
#     33.865,
#     5.09479,
#     1.15879,
#     0.32584,
#     0.102741,
#     0.036,
#     0.75,
#     0.75,
#     0.75,
#     33.865,
#     5.09479,
#     1.15879,
#     0.32584,
#     0.102741,
#     0.036,
#     0.75,
#     0.75,
#     0.75,
# ]

# n = [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1]
# lx = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]
# ly = [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
# lz = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1]

# center = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# total_nprim = 18

# coord = [[0.0, 0.0, 0.0586476414], [0.0, 0.0, 1.4045523587]]

# Norm = {0: phi.NS, 1: phi.NP}
# intLx = np.zeros((total_nprim, total_nprim), dtype=float)
# intLy = np.zeros((total_nprim, total_nprim), dtype=float)
# intLz = np.zeros((total_nprim, total_nprim), dtype=float)
# output = 11

# # !Note: Rgaugeo is put on in the origen
# Rg = [0.0, 0.0, 1.404552358700]
def angmom(coord, gauge, exp, center, lx, ly, lz, output):
    """
    Angular moment integrals, which is a vector

    Args:
        atom (list): list 1d with atoms index
        coord (list): list 2d with coordinates of the atoms
        gauge (list): list 1d with gauge coordinates 
        spatial (list): list with coordinate to evaluate [0:x, 1:y, 2:z]
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
        output (int): Output level for integral calculation

    Return:
        angmom (array): array 2d with atomic integrals
    """

    start = time()
    # Primitive total in the cluster
    total_nprim = len(exp)

    angmom = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    count = 0

    for i in range(total_nprim):

        for j in range(i, total_nprim):

            Px = (
                exp[i] * coord[center[i]][0]
                + exp[j] * coord[center[j]][0]
            )
            Px = Px / (exp[i] + exp[j])

            Py = (
                exp[i] * coord[center[i]][1]
                + exp[j] * coord[center[j]][1]
            )
            Py = Py / (exp[i] + exp[j])

            Pz = (
                exp[i] * coord[center[i]][2]
                + exp[j] * coord[center[j]][2]
            )
            Pz = Pz / (exp[i] + exp[j])

            Xpg = Px - gauge[0]
            Ypg = Py - gauge[1]
            Zpg = Pz - gauge[2]

            sij = E(
                lx[i],
                lx[j],
                0,
                coord[center[i]][0] - coord[center[j]][0],
                exp[i],
                exp[j],
            )

            skl = E(
                ly[i],
                ly[j],
                0,
                coord[center[i]][1] - coord[center[j]][1],
                exp[i],
                exp[j],
            )

            smn = E(
                lz[i],
                lz[j],
                0,
                coord[center[i]][2] - coord[center[j]][2],
                exp[i],
                exp[j],
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
                exp[i],
                exp[j],
            )

            py = 2.0 * exp[j] * E(
                ly[i],
                ly[j] + 1,
                0,
                coord[center[i]][1] - coord[center[j]][1],
                exp[i],
                exp[j],
            ) - ly[j] * E(
                ly[i],
                ly[j] - 1,
                0,
                coord[center[i]][1] - coord[center[j]][1],
                exp[i],
                exp[j],
            )

            zgaugeo = E(
                lz[i],
                lz[j],
                1,
                coord[center[i]][2] - coord[center[j]][2],
                exp[i],
                exp[j],
            )

            pz = 2.0 * exp[j] * E(
                lz[i],
                lz[j] + 1,
                0,
                coord[center[i]][2] - coord[center[j]][2],
                exp[i],
                exp[j],
            ) - lz[j] * E(
                lz[i],
                lz[j] - 1,
                0,
                coord[center[i]][2] - coord[center[j]][2],
                exp[i],
                exp[j],
            )

            intLx[i, j] = (
                Norm[n[i]](exp[i])
                * Norm[n[j]](exp[j])
                * ((ygaugeo + Ypg * skl) * pz - (zgaugeo + Zpg * smn) * py)
                * sij
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )
            intLx[j, i] = -1.0 * intLx[i, j]

            # Real{Ly} = <phi|(zgdx - xgdz)|phi> =
            #  [(zgSmn^0)(dxSij^0) - (xgSij^0)(dzSmn^0)]Skl^0 =
            #  [Smn^1Dij^1 - Sij^1Dmn^1]Skl^0 =
            #  [(E_1^mn + ZpgE0^mn)(2bE_0^ij+1-l2bE0^ij-1) -
            #   (E_1^ij + XpgE0^ij)(2bE_0^mn+1-l2bE0^mn-1)]Skl^0

            xgaugeo = phi.E(
                lx[i],
                lx[j],
                1,
                coord[center[i]][0] - coord[center[j]][0],
                exp[i],
                exp[j],
            )

            px = 2.0 * exp[j] * E(
                lx[i],
                lx[j] + 1,
                0,
                coord[center[i]][0] - coord[center[j]][0],
                exp[i],
                exp[j],
            ) - lx[j] * E(
                lx[i],
                lx[j] - 1,
                0,
                coord[center[i]][0] - coord[center[j]][0],
                exp[i],
                exp[j],
            )

            intLy[i, j] = (
                Norm[n[i]](exp[i])
                * Norm[n[j]](exp[j])
                * ((zgaugeo + Zpg * smn) * px - (xgaugeo + Xpg * sij) * pz)
                * skl
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )
            intLy[j, i] = -1.0 * intLy[i, j]

            # Real{Lz} = <phi|(xgdy - ygdx)|phi> =
            #  [(xgSij^0)(dySkl^0) - (ygSkl^0)(dxSij^0)]Smn^0 =
            #  [Sij^1Dkl^1 - Skl^1Dij^1]Smn^0 =
            #  [(E_1^ij + XpgE0^ij)(2bE_0^kl+1-l2bE0^kl-1) -
            #   (E_1^kl + YpgE0^kl)(2bE_0^ij+1-l2bE0^ij-1)]Smn^0

            intLz[i, j] = (
                Norm[n[i]](exp[i])
                * Norm[n[j]](exp[j])
                * ((xgaugeo + Xpg * sij) * py - (ygaugeo + Ypg * skl) * px)
                * smn
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )
            intLz[j, i] = -1.0 * intLz[i, j]
            # ! Nota: Revisar que tipo de matriz es, summetric, antisymmetric or Cuadrada
            if output > 10 and np.abs(intLz[j, i]) > 1e-2:
                print("int [", j + 1, ",", i + 1, "] : ", intLz[j, i])

    print(" time [s]: ", -start + time.time())
