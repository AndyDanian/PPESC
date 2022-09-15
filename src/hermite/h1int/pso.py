from lib1h import *

def pso(coord, spatial_sym, atom, exp, center, lx, ly, lz, output, dalton_normalization, driver_time):
    """
    Angular moment integrals, which is a vector

    Args:
        coord (list): list 2d with coordinates of the atoms
        spatial_sym (int): spatial symmetry index
        atom (int): atomic index
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
        output (int): Output level for integral calculation
        dalton_normalization (bool): it is used the dalton normalization formule
        drive_time (drv_object): Object to manage the time

    Return:
        pso (array): array 2d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    pso: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    count: int = 0

    """
    Component Selection L = p x r
                        = (zpy-ypz)x + (xpz-zpx)y + (ypx-xpy)z
    where r = r_e - r_k
    """

    der_x_left: int = 0
    der_y_left: int = 0
    der_z_left: int = 0
    der_x_right: int = 0
    der_y_right: int = 0
    der_z_right: int = 0

    r_x_left: int = 0
    r_y_left: int = 0
    r_z_left: int = 0
    r_x_right: int = 0
    r_y_right: int = 0
    r_z_right: int = 0


    if spatial_sym == 0:
        """X Component"""
        der_y_right = 1
        der_z_left = 1
        r_y_left = 1
        r_z_right = 1
        l_right: list = ly
        l_left: list = lz
    elif spatial_sym == 1:
        """Y Componente"""
        der_z_right = 1
        der_x_left = 1
        r_x_right = 1
        r_z_left = 1
        l_right: list = lz
        l_left: list = lx
    elif spatial_sym == 2:
        """Z Componente"""
        der_x_right = 1
        der_y_left = 1
        r_y_right = 1
        r_x_left = 1
        l_right: list = lx
        l_left: list = ly
    else:
        raise ValueError(f"***Error\n\n Component not exist: {spatial_sym}")

    for i in range(total_nprim):

        for j in range(i, total_nprim):

            # pso is a combination of the NELFLD with DPVL
            # (ykdz-zkdy)/rk^3 = yk/rk^3 dz - zk/rk^3 dy =
            # V_ab^010 Dmn^1 - Vab^001 * Dkl^1  (Eq 9.931)

            left_term = 2.0 * exp[j] * nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j] + der_x_left,
                ly[j] + der_y_left,
                lz[j] + der_z_left,
                r_x_left,
                r_y_left,
                r_z_left,
                exp[i],
                exp[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                coord[atom][0],
                coord[atom][1],
                coord[atom][2],
            ) - l_left[j] * nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j] - der_x_left,
                ly[j] - der_y_left,
                lz[j] - der_z_left,
                r_x_left,
                r_y_left,
                r_z_left,
                exp[i],
                exp[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                coord[atom][0],
                coord[atom][1],
                coord[atom][2],
            )

            right_term = 2.0 * exp[j] * nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j] + der_x_right,
                ly[j] + der_y_right,
                lz[j] + der_z_right,
                r_x_right,
                r_y_right,
                r_z_right,
                exp[i],
                exp[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                coord[atom][0],
                coord[atom][1],
                coord[atom][2],
            ) - l_right[j] * nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j] - der_x_right,
                ly[j] - der_y_right,
                lz[j] - der_z_right,
                r_x_right,
                r_y_right,
                r_z_right,
                exp[i],
                exp[j],
                coord[center[i]][0],
                coord[center[i]][1],
                coord[center[i]][2],
                coord[center[j]][0],
                coord[center[j]][1],
                coord[center[j]][2],
                coord[atom][0],
                coord[atom][1],
                coord[atom][2],
            )

            pso[count] = (
                normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
                * (left_term - right_term)
            )
            count += 1
    if output > 10:
        driver_time.add_name_delta_time(name = f"Paramagnetic Spin-Orbit Atomic Integrals, \
        for {spatial_sym} Spatial Symmetry", delta_time  = (time() - start))

    return pso
