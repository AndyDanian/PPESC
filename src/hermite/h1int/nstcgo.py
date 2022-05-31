from lib1h import *

def nstcgo(coord, gauge, spatial_sym, magnetic_component, atom, exp, center, lx, ly, lz, output, dalton_normalization):
    """
    Diamagnetic nuclear shielding tensor atomic integrals

    Args:
        coord (list): list 2d with coordinates of the atoms
        gauge (list): list 1d with gauge coordinates 
        spatial_sym (int): spatial symmetry index
        magnetic_component (int): magnetic component
        atom (int): atomic index
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
        output (int): Output level for integral calculation
        dalton_normalization (bool): it is used the dalton normalization formule

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

    diagonal: bool = False
    if spatial_sym == 0: 
        """X Component"""
        if magnetic_component == 0:
            diagonal = True
            sign: float = 1.0
            r_y_a = 1
            r_y_b = 1
            r_z_c = 1
            coord_ab: int = 1
            coord_c: int = 2
        elif magnetic_component == 1:
            sign: float = -1.0
            coord_ab = 0
            r_x_a = 1
            r_y_b = 1
        elif magnetic_component == 2:
            sign: float = -1.0
            coord_ab = 0
            r_x_a = 1
            r_z_b = 1
    elif spatial_sym == 1:
        """Y Componente"""
        if magnetic_component == 1:
            diagonal = True
            sign: float = 1.0
            r_x_a = 1
            r_x_b = 1
            r_z_c = 1
            coord_ab: int = 0
            coord_c: int = 2
        if magnetic_component == 0:
            coord_ab = 1
            sign: float = -1.0
            r_y_a = 1
            r_x_b = 1
        elif magnetic_component == 2:
            coord_ab = 1
            sign: float = -1.0
            r_y_a = 1
            r_z_b = 1
    elif spatial_sym == 2:
        """Z Componente"""
        if magnetic_component == 2:
            diagonal = True
            sign: float = 1.0
            r_x_a = 1
            r_x_b = 1
            r_y_c = 1
            coord_ab: int = 0
            coord_c: int = 1
        if magnetic_component == 0:
            sign: float = -1.0
            coord_ab = 2
            r_z_a = 1
            r_x_b = 1
        elif magnetic_component == 1:
            sign: float = -1.0
            coord_ab = 2
            r_z_a = 1
            r_y_b = 1
    else:
        raise ValueError(f"***Error\n\n Component not exist: {spatial_sym}")

    nefab: float = 0.0
    nefc: float = 0.0
    for i in range(total_nprim):

        for j in range(i, total_nprim):

            # Nuclear Diamagnetic Shielding
            # <phi|(y-yg)(y-yk)/rk^3 + (z-zg)(z-zk)/rk^3|phi> --> Eqs 9.932, 9.9.18-20
            # *** (y-yg)(y-yk)/r^3_k + (z-zg)(z-zk)/r^3_k ***
            # *** -(x-xg)(y-yk)/r^3_k***
            nefab = sign * (
                nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i],
                    lx[j] + r_x_a,
                    ly[j] + r_y_a,
                    lz[j] + r_z_a,
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
                + (coord[center[j]][coord_ab] - gauge[coord_ab])
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
            )
            if diagonal:
                nefc = (
                    nuclear_attraction(
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
                normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
                * (nefab + nefc)
                * 0.5
            )
            count += 1

    if output > 0:
        print(
            f"\n ***Diamagnetic nuclear shielding tensor atomic integrals \n\
                for {magnetic_component} magnetic component and {spatial_sym} spatial symmetry,\n\
                    time [s]: {time() - start:.6f}"
        )

    return nstcgo
