from libh import *

#################### Calculate the overlap integrals ########################
def overlap(coord, exp, center, lx, ly, lz, output):
    """
    Overlap integrals

    Args:
        charge (list): list 1d of charges
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

    # overlap = [[0 for i in range(total_nprim)] for j in range(total_nprim)]
    overlap = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    count = 0
    for i in range(total_nprim):

        for j in range(i, total_nprim):

            # Build one vector with the elements of a triangular matrix
            # [x x x x x x] = [[x x x][-- x x][-- -- x]]
            # -- are the elements not calculated

            # int -inf inf Hermite Gaussian dx = (pi/p)**1/2 delta_{t0}

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

            overlap[count] = (
                Norm[lx[i] + ly[i] + lz[i]](exp[i])
                * Norm[lx[j] + ly[j] + lz[j]](exp[j])
                * sij
                * skl
                * smn
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )

            count += 1

    if output > 0:
        print(f"\n *** Atomic overlap integrals, time [s]: {time() - start:.6f}")

    return overlap
