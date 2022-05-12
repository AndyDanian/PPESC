from libh import *

############# Calculate the potential one body integrals ########################
def pot(coord, exp, center, lx, ly, lz, output):
    """_summary_

    Potential integrals

    Args:
        coord (list): list 2d with coordinates of the atoms
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
        output (int): Output level for integral calculation

    Return:
        overlap (array): array 2d with atomic integrals
    """

    start = time()
    # Primitive total in the cluster
    total_nprim = len(exp)
    kinetic = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    count = 0
    for i in range(total_nprim):

        for j in range(i, total_nprim):
            # Horizontal Reccurence
            # -0.5*<phi|dxx|phi> = -0.5*(dxxSij^0)Skl^0Smn^0 =
            # -0.5*Dij^2Skl^0Smn^0 =
            # -0.5*(4b^2Sij+2^0 - 2b(2j + 1)Sij^0 + j(j-1)Sij-2^0)Skl^0Smn^0 =
            # -0.5*(4b^2E0^ij+2 - 2b(2j + 1)E0^ij + j(j-1)E0^ij-2)E0^klE0^mn

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

            dxxsij = (
                4.0
                * exp[j]
                * exp[j]
                * (
                    E(
                        lx[i],
                        lx[j] + 2,
                        0,
                        coord[center[i]][0] - coord[center[j]][0],
                        exp[i],
                        exp[j],
                    )
                )
                - 2.0 * exp[j] * (2.0 * lx[j] + 1.0) * sij
                + lx[j]
                * (lx[j] - 1.0)
                * (
                    E(
                        lx[i],
                        lx[j] - 2,
                        0,
                        coord[center[i]][0] - coord[center[j]][0],
                        exp[i],
                        exp[j],
                    )
                )
            )

            dyyskl = (
                4.0
                * exp[j]
                * exp[j]
                * (
                    E(
                        ly[i],
                        ly[j] + 2,
                        0,
                        coord[center[i]][1] - coord[center[j]][1],
                        exp[i],
                        exp[j],
                    )
                )
                - 2.0 * exp[j] * (2.0 * ly[j] + 1.0) * skl
                + ly[j]
                * (ly[j] - 1.0)
                * (
                    E(
                        ly[i],
                        ly[j] - 2,
                        0,
                        coord[center[i]][1] - coord[center[j]][1],
                        exp[i],
                        exp[j],
                    )
                )
            )

            dzzsmn = (
                4.0
                * exp[j]
                * exp[j]
                * (
                    E(
                        lz[i],
                        lz[j] + 2,
                        0,
                        coord[center[i]][2] - coord[center[j]][2],
                        exp[i],
                        exp[j],
                    )
                )
                - 2.0 * exp[j] * (2.0 * lz[j] + 1.0) * smn
                + lz[j]
                * (lz[j] - 1.0)
                * (
                    E(
                        lz[i],
                        lz[j] - 2,
                        0,
                        coord[center[i]][2] - coord[center[j]][2],
                        exp[i],
                        exp[j],
                    )
                )
            )

            kinetic[count] = (
                Norm[lx[i] + ly[i] + lz[i]](exp[i])
                * Norm[lx[j] + ly[j] + lz[j]](exp[j])
                * -0.5
                * (
                    dxxsij * skl * smn
                    + sij * dyyskl * smn
                    + sij * skl * dzzsmn
                )
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )
            count += 1

    if output > 10:
        print("\n *** Atomic kinetic integrals, time [s]: ", time() - start)

    return kinetic, "sym"
