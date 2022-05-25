from libh import *

def psooz(coord, gauge, spatial_sym, magnetic_component, atom, exp, center, lx, ly, lz, output):
    """
    Orbital-Zeeman correction to the paramagnetic spin-orbit atomic integrals

    Agauges:
        coord (list): list 2d with coordinates of the atoms
        gauge (list): list 2d with gauge coordinates
        spatial_sym (int): spatial symmetry index
        magnetic_component (int): magnetic component index
        atom (int): atomic index
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
        output (int): Output level for integral calculation

    Return:
        psooz (array): array 1d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    psooz: list = [0 for i in range(int(total_nprim * total_nprim))]

    count: int = 0
    
    oz_r_x_b: int = 0
    oz_r_y_b: int = 0
    oz_r_z_b: int = 0
    oz_r_x_c: int = 0
    oz_r_y_c: int = 0
    oz_r_z_c: int = 0

    pso_r_x_b: int = 0
    pso_r_y_b: int = 0
    pso_r_z_b: int = 0
    pso_r_x_c: int = 0
    pso_r_y_c: int = 0
    pso_r_z_c: int = 0

    if magnetic_component == 0:
        """X Component"""
        oz_r_y_b = 1
        oz_r_z_c = 1
        coord_b: int = 1
        coord_c: int = 2
        oz_der_l_b: list = ly
        oz_der_l_c: list = lz
    elif magnetic_component == 1:
        """Y Component"""
        oz_r_z_b = 1
        oz_r_x_c = 1
        coord_b: int = 2
        coord_c: int = 0
        oz_der_l_b: list = lz
        oz_der_l_c: list = lx
    elif magnetic_component == 2:
        """Z Component"""
        oz_r_x_b = 1
        oz_r_y_c = 1
        coord_b: int = 0
        coord_c: int = 1
        oz_der_l_b: list = lx
        oz_der_l_c: list = ly

    if spatial_sym == 0: 
        """1 Component"""
        pso_r_y_b = 1
        pso_r_z_c = 1
        pso_der_l_b: list = ly
        pso_der_l_c: list = lz
    elif spatial_sym == 1:
        """2 Componente"""
        pso_r_z_b = 1
        pso_r_x_c = 1
        pso_der_l_b: list = lz
        pso_der_l_c: list = lx
    elif spatial_sym == 2:
        """3 Componente"""
        pso_r_x_b = 1
        pso_r_y_c = 1
        pso_der_l_b: list = lx
        pso_der_l_c: list = ly
    else:
        raise ValueError(f"***Error\n\n Component not exist: {spatial_sym}")

    for i in range(total_nprim):

        for j in range(total_nprim):

            # (ygdz-zgdy)(ykdz-zkdy)/rk^3 =
            # (ygdz) yk/rk^3 dz - (ygdz) zk/rk^3 dy - (zgdy) yk/rk^3 dz + (zgdy) zk/rk^3 dy=
            # (yk + ky - yg) Dmn^1 V_ab^010 Dmn^1 - (yk + ky - yg) Dmn^1 Vab^001 Dkl^1
            # (zk + kz - zg) Dkl^1 V_ab^010 Dmn^1 - (zk + kz - zg) Dkl^1 Vab^001 Dkl^1 (Eq 9.931)
            # (zgdx-xgdz)(ykdz-zkdy)/rk^3

            igdj_jgdi_idj_jdi: float = 0.0
            for k in range(8):
                if k == 0:
                    l_x_i_a, l_y_i_a, l_z_i_a = oz_r_x_b + oz_r_x_c, oz_r_y_b + oz_r_y_c, oz_r_z_b + oz_r_z_c
                    l_x_i_b, l_y_i_b, l_z_i_b = oz_r_x_b - oz_r_x_c, oz_r_y_b - oz_r_y_c, oz_r_z_b - oz_r_z_c
                    l_x_j, l_y_j, l_z_j = pso_r_x_c, pso_r_y_c, pso_r_z_c
                    r_x, r_y, r_z = pso_r_x_b, pso_r_y_b, pso_r_z_b
                    coef: float = 1.0
                    der_l_left: float = float(oz_der_l_c[i])
                    der_l_right: float = float(pso_der_l_c[j])
                elif k == 1:
                    l_x_i_a, l_y_i_a, l_z_i_a = oz_r_x_c, oz_r_y_c, oz_r_z_c
                    l_x_i_b, l_y_i_b, l_z_i_b = -oz_r_x_c, -oz_r_y_c, -oz_r_z_c
                    l_x_j, l_y_j, l_z_j = pso_r_x_c, pso_r_y_c, pso_r_z_c
                    r_x, r_y, r_z = pso_r_x_b, pso_r_y_b, pso_r_z_b
                    coef: float = 1.0 * (coord[center[i]][coord_b] - gauge[coord_b])
                    der_l_left: float = float(oz_der_l_c[i])
                    der_l_right: float = float(pso_der_l_c[j])
                elif k == 2:
                    l_x_i_a, l_y_i_a, l_z_i_a = oz_r_x_b + oz_r_x_c, oz_r_y_b + oz_r_y_c, oz_r_z_b + oz_r_z_c
                    l_x_i_b, l_y_i_b, l_z_i_b = -oz_r_x_b + oz_r_x_c, -oz_r_y_b + oz_r_y_c, -oz_r_z_b + oz_r_z_c
                    l_x_j, l_y_j, l_z_j = pso_r_x_c, pso_r_y_c, pso_r_z_c
                    r_x, r_y, r_z = pso_r_x_b, pso_r_y_b, pso_r_z_b
                    coef: float = -1.0
                    der_l_left: float = float(oz_der_l_b[i])
                    der_l_right: float = float(pso_der_l_c[j])
                elif k == 3:
                    l_x_i_a, l_y_i_a, l_z_i_a = oz_r_x_b, oz_r_y_b, oz_r_z_b
                    l_x_i_b, l_y_i_b, l_z_i_b = -oz_r_x_b, -oz_r_y_b, -oz_r_z_b
                    l_x_j, l_y_j, l_z_j = pso_r_x_c, pso_r_y_c, pso_r_z_c
                    r_x, r_y, r_z = pso_r_x_b, pso_r_y_b, pso_r_z_b
                    coef: float = -1.0 * (coord[center[i]][coord_c] - gauge[coord_c])
                    der_l_left: float = float(oz_der_l_b[i])
                    der_l_right: float = float(pso_der_l_c[j])
                elif k == 4:
                    l_x_i_a, l_y_i_a, l_z_i_a = oz_r_x_b + oz_r_x_c, oz_r_y_b + oz_r_y_c, oz_r_z_b + oz_r_z_c
                    l_x_i_b, l_y_i_b, l_z_i_b = oz_r_x_b - oz_r_x_c, oz_r_y_b - oz_r_y_c, oz_r_z_b - oz_r_z_c
                    l_x_j, l_y_j, l_z_j = pso_r_x_b, pso_r_y_b, pso_r_z_b
                    r_x, r_y, r_z = pso_r_x_c, pso_r_y_c, pso_r_z_c
                    coef: float = -1.0
                    der_l_left: float = float(oz_der_l_c[i])
                    der_l_right: float = float(pso_der_l_b[j])
                elif k == 5:
                    l_x_i_a, l_y_i_a, l_z_i_a = oz_r_x_c, oz_r_y_c, oz_r_z_c
                    l_x_i_b, l_y_i_b, l_z_i_b = -oz_r_x_c, -oz_r_y_c, -oz_r_z_c
                    l_x_j, l_y_j, l_z_j = pso_r_x_b, pso_r_y_b, pso_r_z_b
                    r_x, r_y, r_z = pso_r_x_c, pso_r_y_c, pso_r_z_c
                    coef: float = -1.0 * (coord[center[i]][coord_b] - gauge[coord_b])
                    der_l_left: float = float(oz_der_l_c[i])
                    der_l_right: float = float(pso_der_l_b[j])
                elif k == 6:
                    l_x_i_a, l_y_i_a, l_z_i_a = oz_r_x_b + oz_r_x_c, oz_r_y_b + oz_r_y_c, oz_r_z_b + oz_r_z_c
                    l_x_i_b, l_y_i_b, l_z_i_b = -oz_r_x_b + oz_r_x_c, -oz_r_y_b + oz_r_y_c, -oz_r_z_b + oz_r_z_c
                    l_x_j, l_y_j, l_z_j = pso_r_x_b, pso_r_y_b, pso_r_z_b
                    r_x, r_y, r_z = pso_r_x_c, pso_r_y_c, pso_r_z_c
                    coef: float = 1.0
                    der_l_left: float = float(oz_der_l_b[i])
                    der_l_right: float = float(pso_der_l_b[j])
                elif k == 7:
                    l_x_i_a, l_y_i_a, l_z_i_a = oz_r_x_b, oz_r_y_b, oz_r_z_b
                    l_x_i_b, l_y_i_b, l_z_i_b = -oz_r_x_b, -oz_r_y_b, -oz_r_z_b
                    l_x_j, l_y_j, l_z_j = pso_r_x_b, pso_r_y_b, pso_r_z_b
                    r_x, r_y, r_z = pso_r_x_c, pso_r_y_c, pso_r_z_c
                    coef: float = 1.0 * (coord[center[i]][coord_c] - gauge[coord_c])
                    der_l_left: float = float(oz_der_l_b[i])
                    der_l_right: float = float(pso_der_l_b[j])

                igdj_jgdi_idj_jdi += coef * (
                    4.0
                    * exp[i]
                    * exp[j]
                    * nuclear_attraction(
                        lx[i] + l_x_i_a,
                        ly[i] + l_y_i_a,
                        lz[i] + l_z_i_a,
                        lx[j] + l_x_j,
                        ly[j] + l_y_j,
                        lz[j] + l_z_j,
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
                    - 2.0
                    * exp[i]
                    * der_l_right
                    * nuclear_attraction(
                        lx[i] + l_x_i_a,
                        ly[i] + l_y_i_a,
                        lz[i] + l_z_i_a,
                        lx[j] - l_x_j,
                        ly[j] - l_y_j,
                        lz[j] - l_z_j,
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
                    - 2.0
                    * exp[j]
                    * der_l_left
                    * nuclear_attraction(
                        lx[i] + l_x_i_b,
                        ly[i] + l_y_i_b,
                        lz[i] + l_z_i_b,
                        lx[j] + l_x_j,
                        ly[j] + l_y_j,
                        lz[j] + l_z_j,
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
                    + der_l_left
                    * der_l_right
                    * nuclear_attraction(
                        lx[i] + l_x_i_b,
                        ly[i] + l_y_i_b,
                        lz[i] + l_z_i_b,
                        lx[j] - l_x_j,
                        ly[j] - l_y_j,
                        lz[j] - l_z_j,
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

            psooz[count] = (
                normalization(lx[i], ly[i], lz[i], exp[i])
                * normalization(lx[j], ly[j], lz[j], exp[j])
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
                * igdj_jgdi_idj_jdi
            )
            count += 1

    if output > 10:
        print(f"\n ***Orbital-Zeeman correction to the paramagnetic spin-orbit atomic integrals,\n\
        for {magnetic_component} magnetic component, {spatial_sym} spatial symmetry, and {atom} atom,\n\
        time [s]: {time() - start:.6f}")

    return psooz