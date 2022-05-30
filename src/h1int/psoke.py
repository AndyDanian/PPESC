from libh import *
def psoke(coord, spatial_sym, atom, exp, center, lx, ly, lz, output, dalton_normalization):
    """
    Kinetic energy correction to the paramagnetic spin-orbit atomic integrals

    Agauges:
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

    Return:
        psoke (array): array 1d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    psoke: list = [0 for i in range(int(total_nprim * total_nprim))]

    count: int = 0
    
    r_x_b: int = 0
    r_y_b: int = 0
    r_z_b: int = 0
    r_x_c: int = 0
    r_y_c: int = 0
    r_z_c: int = 0

    if spatial_sym == 0: 
        """X Component"""
        r_y_b = 1
        r_z_c = 1
        der_l_b: list = ly
        der_l_c: list = lz
    elif spatial_sym == 1:
        """Y Componente"""
        r_z_b = 1
        r_x_c = 1
        der_l_b: list = lz
        der_l_c: list = lx
    elif spatial_sym == 2:
        """Z Componente"""
        r_x_b = 1
        r_y_c = 1
        der_l_b: list = lx
        der_l_c: list = ly
    else:
        raise ValueError(f"***Error\n\n Component not exist: {spatial_sym}")

    lap: list = [(2,0,0),(0,2,0),(0,0,2)]
    for i in range(total_nprim):

        for j in range(total_nprim):

            idj_jdi: list = []
            for k in range(2):
                if k == 0:
                    l_x, l_y, l_z = r_x_c, r_y_c, r_z_c
                    r_x, r_y, r_z = r_x_b, r_y_b, r_z_b
                    sign: float = 1.0
                    der_l: float = float(der_l_c[j])
                else:
                    l_x, l_y, l_z = r_x_b, r_y_b, r_z_b
                    r_x, r_y, r_z = r_x_c, r_y_c, r_z_c
                    sign: float = -1.0
                    der_l: float = float(der_l_b[j])
                idj_jdi.append(
                    2.0 * exp[j] * nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i],
                    lx[j] + l_x,
                    ly[j] + l_y,
                    lz[j] + l_z,
                    r_x,
                    r_y,
                    r_z,
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
                ) - 
                der_l * nuclear_attraction(
                    lx[i],
                    ly[i],
                    lz[i],
                    lx[j] - l_x,
                    ly[j] - l_y,
                    lz[j] - l_z,
                    r_x,
                    r_y,
                    r_z,
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
                ))
            
            # Todo bien hasta áca
            lap_idj_jdi: float = 0.0
            idj_jdi_lap: float = 0.0
            for d2x, d2y, d2z in lap:
                if d2x == 2:
                    lap_l_i: float = float(lx[i])
                    lap_l_j: float = float(lx[j])
                elif d2y == 2:
                    lap_l_i: float = float(ly[i])
                    lap_l_j: float = float(ly[j])
                elif d2z == 2:
                    lap_l_i: float = float(lz[i])
                    lap_l_j: float = float(lz[j])

                for k in range(2):
                    if k == 0:
                        l_x, l_y, l_z = r_x_c, r_y_c, r_z_c
                        r_x, r_y, r_z = r_x_b, r_y_b, r_z_b
                        sign: float = 1.0
                        der_l: float = float(der_l_c[j])
                    else:
                        l_x, l_y, l_z = r_x_b, r_y_b, r_z_b
                        r_x, r_y, r_z = r_x_c, r_y_c, r_z_c
                        sign: float = -1.0
                        der_l: float = float(der_l_b[j])
                # ! Terms of dxxLxyz/rk^3
                    lap_idj_jdi += sign*(
                    4.0
                    * exp[i]
                    * exp[i]
                    * (
                        2.0
                        * exp[j]
                        * nuclear_attraction(
                            lx[i] + d2x,
                            ly[i] + d2y,
                            lz[i] + d2z,
                            lx[j] + l_x,
                            ly[j] + l_y,
                            lz[j] + l_z,
                            r_x,
                            r_y,
                            r_z,
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
                        - der_l
                        * nuclear_attraction(
                            lx[i] + d2x,
                            ly[i] + d2y,
                            lz[i] + d2z,
                            lx[j] - l_x,
                            ly[j] - l_y,
                            lz[j] - l_z,
                            r_x,
                            r_y,
                            r_z,
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
                    - 2.0 * exp[i] * (2.0 * lap_l_i + 1.0) * idj_jdi[k]
                    + lap_l_i
                    * (lap_l_i - 1.0)
                    * (
                        2.0
                        * exp[j]
                        * nuclear_attraction(
                            lx[i] - d2x,
                            ly[i] - d2y,
                            lz[i] - d2z,
                            lx[j] + l_x,
                            ly[j] + l_y,
                            lz[j] + l_z,
                            r_x,
                            r_y,
                            r_z,
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
                        - der_l
                        * nuclear_attraction(
                            lx[i] - d2x,
                            ly[i] - d2y,
                            lz[i] - d2z,
                            lx[j] - r_x,
                            ly[j] - r_y,
                            lz[j] - r_z,
                            r_x,
                            r_y,
                            r_z,
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
                        ) ))
            # * nabla Real{Lx} + Real{Lx} nabla

            psoke[count] = (
                normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * (lap_idj_jdi)
                * 0.5
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
            )
            count += 1
    if output > 10:
        print(f"\n ***Kinetic energy correction to the paramagnetic spin-orbit atomic integrals,\n\
        for {spatial_sym} spatial symmetry, time [s]: {time() - start:.6f}")

    return psoke
