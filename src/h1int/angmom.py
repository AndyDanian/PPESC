from libh import *

def angmom(coord, gauge, magnetic_component, exp, center, lx, ly, lz, output):
    """
    Angular moment integrals, which is a vector

    Args:
        coord (list): list 2d with coordinates of the atoms
        gauge (list): list 1d with gauge coordinates 
        magnetic_component (int): magnetic component
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
        output (int): Output level for integral calculation

    Return:
        angmom (array): array 2d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    angmom: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    count: int = 0

    """
    Component Selection L = p x r 
                          = (zpy-ypz)x + (xpz-zpx)y + (ypx-xpy)z
    where r = r_e - r_gauge
    """

    if magnetic_component == 0: 
        """X Component"""
        left_coord: int = 1
        right_coord: int = 2
        spatial_l: list = lx
        left_l: list = ly
        right_l: list = lz
    elif magnetic_component == 1:
        """Y Componente"""
        left_coord: int = 2
        right_coord: int = 0
        spatial_l: list = ly
        left_l: list = lz
        right_l: list = lx
    elif magnetic_component == 2:
        """Z Componente"""
        left_coord: int = 0
        right_coord: int = 1
        spatial_l: list = lz
        left_l: list = lx
        right_l: list = ly
    else:
        raise ValueError(f"***Error\n\n Component not exist: {magnetic_component}")

    for i in range(total_nprim):

        for j in range(i, total_nprim):

            # Gaussian Center
            left_Pxyz: float = (
                exp[i] * coord[center[i]][left_coord]
                + exp[j] * coord[center[j]][left_coord]
            )
            left_Pxyz = left_Pxyz / (exp[i] + exp[j])
            right_Pxyz: float = (
                exp[i] * coord[center[i]][right_coord]
                + exp[j] * coord[center[j]][right_coord]
            )
            right_Pxyz = right_Pxyz / (exp[i] + exp[j])

            left_rg: float = left_Pxyz - gauge[left_coord]
            right_rg: float = right_Pxyz - gauge[right_coord]
    
            spatial_s: float = E(
                spatial_l[i],
                spatial_l[j],
                0,
                coord[center[i]][magnetic_component] - coord[center[j]][magnetic_component],
                exp[i],
                exp[j],
            )

            left_s: float = E(
                left_l[i],
                left_l[j],
                0,
                coord[center[i]][left_coord] - coord[center[j]][left_coord],
                exp[i],
                exp[j],
            )

            right_s: float = E( #z
                right_l[i],
                right_l[j],
                0,
                coord[center[i]][right_coord] - coord[center[j]][right_coord],
                exp[i],
                exp[j],
            )

            # Real{Lx} = <phi|(ygdz - zgdy)|phi> =
            #  [(ygSkl^0)(dzSmn^0) - (zgSmn^0)(dySkl^0)]Sij^0 =
            #  [Skl^1Dmn^1 - Smn^1Dkl^1]Sij^0 =
            #  [(E_1^kl + YpgE0^kl)(2bE_0^mn+1-l2bE0^mn-1) -
            #   (E_1^mn + ZpgE0^mn)(2bE_0^kl+1-l2bE0^kl-1)]Sij^0
            left_r: float = E(
                left_l[i],
                left_l[j],
                1,
                coord[center[i]][left_coord] - coord[center[j]][left_coord],
                exp[i],
                exp[j],
            )

            right_r: float = E(
                right_l[i],
                right_l[j],
                1,
                coord[center[i]][right_coord] - coord[center[j]][right_coord],
                exp[i],
                exp[j],
            )

            right_p: float = 2.0 * exp[j] * E(
                left_l[i],
                left_l[j] + 1,
                0,
                coord[center[i]][left_coord] - coord[center[j]][left_coord],
                exp[i],
                exp[j],
            ) - left_l[j] * E(
                left_l[i],
                left_l[j] - 1,
                0,
                coord[center[i]][left_coord] - coord[center[j]][left_coord],
                exp[i],
                exp[j],
            )

            left_p: float = 2.0 * exp[j] * E(
                right_l[i],
                right_l[j] + 1,
                0,
                coord[center[i]][right_coord] - coord[center[j]][right_coord],
                exp[i],
                exp[j],
            ) - right_l[j] * E(
                right_l[i],
                right_l[j] - 1,
                0,
                coord[center[i]][right_coord] - coord[center[j]][right_coord],
                exp[i],
                exp[j],
            )
    
            angmom[count] = (
                -Norm[lx[i] + ly[i] + lz[i]](exp[i])
                * Norm[lx[j] + ly[j] + lz[j]](exp[j])
                * ((left_r + left_rg * left_s) * left_p -
                   (right_r + right_rg * right_s) * right_p)
                * spatial_s
                * np.power(np.pi / (exp[i] + exp[j]), 1.5)
            )
            count += 1

    if output > 10:
        print(f"\n *** Angular momentum atomic integrals,\
        component {magnetic_component}, time [s]: {time() - start:.6f}")

    return angmom