from libh import *

def angmom(coord, gauge, spatial_sym, exp, center, lx, ly, lz, output):
    """
    Angular moment integrals, which is a vector

    Args:
        coord (list): list 2d with coordinates of the atoms
        gauge (list): list 1d with gauge coordinates 
        spatial_sym (list): list with coordinate to evaluate [0:x, 1:y, 2:z]
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

    count: int = 0

    """
    Component Selection L = p x r 
                          = (zpy-ypz)x + (xpz-zpx)y + (ypx-xpy)z
    """

    if spatial_sym == 0: 
        """X Component"""
        left_coord: int = 1
        right_coord: int = 2
        left_l: list = ly
        right_l: list = lz
    elif spatial_sym == 1:
        """Y Componente"""
        left_coord: int = 2
        right_coord: int = 0
        left_l: list = lz
        right_l: list = lx
    elif spatial_sym == 2:
        """Z Componente"""
        left_coord: int = 0
        right_coord: int = 1
        left_l: list = lx
        right_l: list = ly
    else:
        raise ValueError(f"***Error\n\n Component not exist: {spatial_sym}")


    for i in range(total_nprim):

        for j in range(i, total_nprim):

            # Px = (
            #     exp[i] * coord[center[i]][0]
            #     + exp[j] * coord[center[j]][0]
            # )
            # Px = Px / (exp[i] + exp[j])

            # Py = (
            #     exp[i] * coord[center[i]][1]
            #     + exp[j] * coord[center[j]][1]
            # )
            # Py = Py / (exp[i] + exp[j])

            # Pz = (
            #     exp[i] * coord[center[i]][2]
            #     + exp[j] * coord[center[j]][2]
            # )
            # Pz = Pz / (exp[i] + exp[j])

            Pxyz = (
                exp[i] * coord[center[i]][spatial_sym]
                + exp[j] * coord[center[j]][spatial_sym]
            )
            Pxyz = Pxyz / (exp[i] + exp[j])

            # Xpg = Px - gauge[0]
            # Ypg = Py - gauge[1]
            # Zpg = Pz - gauge[2]
            rg = Pxyz - gauge[spatial_sym]

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
            # ygaugeo = E(
            #     ly[i],
            #     ly[j],
            #     1,
            #     coord[center[i]][1] - coord[center[j]][1],
            #     exp[i],
            #     exp[j],
            # )

            # xgaugeo = E(
            #     lx[i],
            #     lx[j],
            #     1,
            #     coord[center[i]][0] - coord[center[j]][0],
            #     exp[i],
            #     exp[j],
            # )

            # zgaugeo = E(
            #     lz[i],
            #     lz[j],
            #     1,
            #     coord[center[i]][2] - coord[center[j]][2],
            #     exp[i],
            #     exp[j],
            # )


            left_r = E(
                left_l[i],
                left_l[j],
                1,
                coord[center[i]][left_coord] - coord[center[j]][left_coord],
                exp[i],
                exp[j],
            )

            right_r = E(
                right_l[i],
                right_l[j],
                1,
                coord[center[i]][right_coord] - coord[center[j]][right_coord],
                exp[i],
                exp[j],
            )


            # py = 2.0 * exp[j] * E(
            #     ly[i],
            #     ly[j] + 1,
            #     0,
            #     coord[center[i]][1] - coord[center[j]][1],
            #     exp[i],
            #     exp[j],
            # ) - ly[j] * E(
            #     ly[i],
            #     ly[j] - 1,
            #     0,
            #     coord[center[i]][1] - coord[center[j]][1],
            #     exp[i],
            #     exp[j],
            # )


            # pz = 2.0 * exp[j] * E(
            #     lz[i],
            #     lz[j] + 1,
            #     0,
            #     coord[center[i]][2] - coord[center[j]][2],
            #     exp[i],
            #     exp[j],
            # ) - lz[j] * E(
            #     lz[i],
            #     lz[j] - 1,
            #     0,
            #     coord[center[i]][2] - coord[center[j]][2],
            #     exp[i],
            #     exp[j],
            # )

            right_p = 2.0 * exp[j] * E(
                right_l[i],
                right_l[j] + 1,
                0,
                coord[center[i]][left_coord] - coord[center[j]][left_coord],
                exp[i],
                exp[j],
            ) - right_l[j] * E(
                right_l[i],
                right_l[j] - 1,
                0,
                coord[center[i]][left_coord] - coord[center[j]][left_coord],
                exp[i],
                exp[j],
            )

            left_p = 2.0 * exp[j] * E(
                left_l[i],
                left_l[j] + 1,
                0,
                coord[center[i]][right_coord] - coord[center[j]][right_coord],
                exp[i],
                exp[j],
            ) - left_l[j] * E(
                left_l[i],
                left_l[j] - 1,
                0,
                coord[center[i]][2] - coord[center[right_coord]][right_coord],
                exp[i],
                exp[j],
            )


            angmom[count] = (
                Norm[lx[i] + ly[i] + lz[i]](exp[i])
                * Norm[lx[j] + ly[j] + lz[j]](exp[j])
                * ((left_r + rg * skl) * left_p - (right_r + rg * smn) * right_p)
                * sij
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )
            count += 1

    if output > 10:
        print(f"\n *** Atomic angular momentum integrals time [s]: {time() - start:.6f}")

    return angmom