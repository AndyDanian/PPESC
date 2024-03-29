from lib1h import *

#################### Calculate the overlap integrals ########################
def overlap(coord, exp, center, lx, ly, lz, output, dalton_normalization, driver_time):
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
        dalton_normalization (bool): it is used the dalton normalization formule
        drive_time (drv_object): Object to manage the time

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

            sij = hermite_coefficient(
                lx[i],
                lx[j],
                0,
                coord[center[i]][0] - coord[center[j]][0],
                exp[i],
                exp[j],
            )

            skl = hermite_coefficient(
                ly[i],
                ly[j],
                0,
                coord[center[i]][1] - coord[center[j]][1],
                exp[i],
                exp[j],
            )

            smn = hermite_coefficient(
                lz[i],
                lz[j],
                0,
                coord[center[i]][2] - coord[center[j]][2],
                exp[i],
                exp[j],
            )

            overlap[count] = (
                normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * sij
                * skl
                * smn
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )
            count += 1

    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Overlap Atomic Integrals", delta_time=(time() - start)
        )

    return overlap


if __name__ == "__main__":
    # 6-311++G**
    s = overlap(
        coord=[[0.0, 0.0, 0.0586476414], [0.0, 0.0, 1.4045523587]],
        exp=[
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
        ],
        center=[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        lx=[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        ly=[0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        lz=[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        output=11,
        dalton_normalization=False,
    )

    print("Overlap : ", s)
