from libh import *

def nstcgo(coord, gauge, spatial_sym, magnetic_component, atom, exp, center, lx, ly, lz, output):
    """
    Diamagnetic nuclear shielding tensor atomic integrals

    Args:
        coord (list): list 2d with coordinates of the atoms
        gauge (list): list 1d with gauge coordinates 
        magnetic_component (int): magnetic component
        atom (int): atomic index
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

    nstcgo: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    count: int = 0

    l_x: int = 0
    l_y: int = 0
    l_z: int = 0
    r_x_a: int = 0
    r_y_a: int = 0
    r_z_a: int = 0
    r_x_b: int = 0
    r_y_b: int = 0
    r_z_b: int = 0
    r_x_c: int = 0
    r_y_c: int = 0
    r_z_c: int = 0

    out_diagonal: bool = False
    if spatial_sym == 0: 
        """X Component"""
        r_y_b = 1
        r_z_c = 1
        coord_b: int = 1
        coord_c: int = 2
        if magnetic_component == 1:
            out_diagonal = True
            coord_a = 0
            l_x = 1
            r_y_a = 1
        elif magnetic_component == 2:
            out_diagonal = True
            coord_a = 0
            l_x = 1
            r_z_a = 1
    elif spatial_sym == 1:
        """Y Componente"""
        r_x_b = 1
        r_z_c = 1
        coord_b: int = 0
        coord_c: int = 2
        if magnetic_component == 0:
            out_diagonal = True
            coord_a = 1
            l_y = 1
            r_x_a = 1
        elif magnetic_component == 2:
            out_diagonal = True
            coord_a = 1
            l_y = 1
            r_z_a = 1
    elif spatial_sym == 2:
        """Z Componente"""
        r_x_b = 1
        r_y_c = 1
        coord_b: int = 0
        coord_c: int = 1
        if magnetic_component == 0:
            out_diagonal = True
            coord_a = 2
            l_z = 1
            r_x_a = 1
        elif magnetic_component == 1:
            out_diagonal = True
            coord_a = 2
            l_z = 1
            r_y_a = 1
    else:
        raise ValueError(f"***Error\n\n Component not exist: {spatial_sym}")

    nefa: float = 0.0
    nefbc: float = 0.0
    for i in range(total_nprim):

        for j in range(i, total_nprim):

            # Nuclear Diamagnetic Shielding
            # <phi|(y-yg)(y-yk)/rk^3 + (z-zg)(z-zk)/rk^3|phi> --> Eqs 9.932, 9.9.18-20
            # *** (y-yg)(y-yk)/r^3_k + (z-zg)(z-zk)/r^3_k ***
            # *** -(x-xg)(y-yk)/r^3_k***
            if out_diagonal:
                nefa = -1.0 * (
                nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i],
                    lx[j] + l_x,
                    ly[j] + l_y,
                    lz[j] + l_z,
                    r_x_a,
                    r_y_a,
                    r_z_a,
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
                + (coord[center[j]][coord_a] - gauge[coord_a])
                * nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i],
                    lx[j],
                    ly[j],
                    lz[j],
                    r_x_a,
                    r_y_a,
                    r_z_a,
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
            )
            else:
                # Diagonals terms
                nefbc = (
                    nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + r_x_b,
                        ly[j] + r_y_b,
                        lz[j] + r_z_b,
                        r_x_b,
                        r_y_b,
                        r_z_b,
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
                    + (coord[center[j]][coord_b] - gauge[coord_b])
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        r_x_b,
                        r_y_b,
                        r_z_b,
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
                    + nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j] + r_x_c,
                        ly[j] + r_y_c,
                        lz[j] + r_z_c,
                        r_x_c,
                        r_y_c,
                        r_z_c,
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
                    + (coord[center[j]][coord_c] - gauge[coord_c])
                    * nuclear_attraction(
                        lx[i],
                        ly[i],
                        lz[i],
                        lx[j],
                        ly[j],
                        lz[j],
                        r_x_c,
                        r_y_c,
                        r_z_c,
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
                )
            # ! Falta poner el -1^e+g+f en la recurrencia de nuclear_attraction
            nstcgo[count] = (
                Norm[lx[i] + ly[i] + lz[i]](exp[i])
                * Norm[lx[j] + ly[j] + lz[j]](exp[j])
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
                * (nefa + nefbc)
                * 0.5
            )
            count += 1

    if output > 0:
        print(
            f"\n ***Diamagnetic nuclear shielding tensor\
                 atomic integrals for {magnetic_component}, time [s]: {time() - start:.6f}"
        )

    return nstcgo
