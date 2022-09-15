from lib1h import *

def pnstcgop(coord, gauge, spatial_sym, magnetic_component, atom,
            exp, center, lx, ly, lz, output, dalton_normalization, driver_time):
    """
    Gradient of diamagnetic nuclear shielding tensor atomic integrals

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
        drive_time (drv_object): Object to manage the time

    Return:
        pnstcgop (array): array 1d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    pnstcgop: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    count: int = 0

    # p
    ml_i_a: int = 0
    ml_k_a: int = 0
    ml_m_a: int = 0
    ml_i_b: int = 0
    ml_k_b: int = 0
    ml_m_b: int = 0
    ml_a: list = []
    ml_b: list = []
    # A²
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
            ml_i_a = 1
            ml_i_b = 1
            r_y_a = 1
            r_y_b = 1
            r_z_c = 1
            ml_a = lx
            ml_b = lx
            coord_ab: int = 1
            coord_c: int = 2
        elif magnetic_component == 1:
            sign: float = -1.0
            ml_k_a = 1
            ml_i_b = 1
            ml_a = ly
            ml_b = lx
            coord_ab = 0
            r_x_a = 1
            r_y_b = 1
        elif magnetic_component == 2:
            sign: float = -1.0
            ml_m_a = 1
            ml_i_b = 1
            ml_a = lz
            ml_b = lx
            coord_ab = 0
            r_x_a = 1
            r_z_b = 1
    elif spatial_sym == 1:
        """Y Componente"""
        if magnetic_component == 1:
            diagonal = True
            sign: float = 1.0
            ml_i_a = 1
            ml_k_b = 1
            r_x_a = 1
            r_x_b = 1
            r_z_c = 1
            ml_a = lx
            ml_b = ly
            coord_ab: int = 0
            coord_c: int = 2
        if magnetic_component == 0:
            coord_ab = 1
            sign: float = -1.0
            ml_k_a = 1
            ml_k_b = 1
            r_y_a = 1
            r_x_b = 1
            ml_a = ly
            ml_b = ly
        elif magnetic_component == 2:
            coord_ab = 1
            sign: float = -1.0
            ml_m_a = 1
            ml_k_b = 1
            r_y_a = 1
            r_z_b = 1
            ml_a = lz
            ml_b = ly
    elif spatial_sym == 2:
        """Z Componente"""
        if magnetic_component == 2:
            diagonal = True
            sign: float = 1.0
            ml_i_a = 1
            ml_m_b = 1
            r_x_a = 1
            r_x_b = 1
            r_y_c = 1
            coord_ab: int = 0
            coord_c: int = 1
            ml_a = lx
            ml_b = lz
        if magnetic_component == 0:
            sign: float = -1.0
            ml_k_a = 1
            ml_m_b = 1
            coord_ab = 2
            r_z_a = 1
            r_x_b = 1
            ml_a = ly
            ml_b = lz
        elif magnetic_component == 1:
            sign: float = -1.0
            ml_m_a = 1
            ml_m_b = 1
            coord_ab = 2
            r_z_a = 1
            r_y_b = 1
            ml_a = lz
            ml_b = lz
    else:
        raise ValueError(f"***Error\n\n Component not exist: {spatial_sym}")

    nefab: float = 0.0
    nefc: float = 0.0
    for i in range(total_nprim):

        for j in range(i, total_nprim):

            # (di psi) A²_ij (dj psi)
            nefab = sign * (
                4.0*exp[i]*exp[j]*
                (nuclear_attraction(
                    lx[i] + ml_i_a,
                    ly[i] + ml_k_a,
                    lz[i] + ml_m_a,
                    lx[j] + r_x_a + ml_i_b,
                    ly[j] + r_y_a + ml_k_b,
                    lz[j] + r_z_a + ml_m_b,
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
                    lx[i] + ml_i_a,
                    ly[i] + ml_k_a,
                    lz[i] + ml_m_a,
                    lx[j] + ml_i_b,
                    ly[j] + ml_k_b,
                    lz[j] + ml_m_b,
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
                - 2.0*exp[i]*ml_b[j]*(
                    nuclear_attraction(
                    lx[i] + ml_i_a,
                    ly[i] + ml_k_a,
                    lz[i] + ml_m_a,
                    lx[j] + r_x_a - ml_i_b,
                    ly[j] + r_y_a - ml_k_b,
                    lz[j] + r_z_a - ml_m_b,
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
                    lx[i] + ml_i_a,
                    ly[i] + ml_k_a,
                    lz[i] + ml_m_a,
                    lx[j] - ml_i_b,
                    ly[j] - ml_k_b,
                    lz[j] - ml_m_b,
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
                - 2.0*exp[j]*ml_a[i]*(
                    nuclear_attraction(
                    lx[i] - ml_i_a,
                    ly[i] - ml_k_a,
                    lz[i] - ml_m_a,
                    lx[j] + r_x_a + ml_i_b,
                    ly[j] + r_y_a + ml_k_b,
                    lz[j] + r_z_a + ml_m_b,
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
                    lx[i] - ml_i_a,
                    ly[i] - ml_k_a,
                    lz[i] - ml_m_a,
                    lx[j] + ml_i_b,
                    ly[j] + ml_k_b,
                    lz[j] + ml_m_b,
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
                + ml_a[i]*ml_b[j]*(
                    nuclear_attraction(
                    lx[i] - ml_i_a,
                    ly[i] - ml_k_a,
                    lz[i] - ml_m_a,
                    lx[j] + r_x_a - ml_i_b,
                    ly[j] + r_y_a - ml_k_b,
                    lz[j] + r_z_a - ml_m_b,
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
                    lx[i] - ml_i_a,
                    ly[i] - ml_k_a,
                    lz[i] - ml_m_a,
                    lx[j] - ml_i_b,
                    ly[j] - ml_k_b,
                    lz[j] - ml_m_b,
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
            )
            if diagonal:
                nefc = (
                    4.0*exp[i]*exp[j]*(nuclear_attraction(
                        lx[i] + ml_i_a,
                        ly[i] + ml_k_a,
                        lz[i] + ml_m_a,
                        lx[j] + r_x_c + ml_i_b,
                        ly[j] + r_y_c + ml_k_b,
                        lz[j] + r_z_c + ml_m_b,
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
                        lx[i] + ml_i_a,
                        ly[i] + ml_k_a,
                        lz[i] + ml_m_a,
                        lx[j] + ml_i_b,
                        ly[j] + ml_k_b,
                        lz[j] + ml_m_b,
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
                    ))
                        - 2.0*exp[i]*ml_b[j]*(
                            nuclear_attraction(
                            lx[i] + ml_i_a,
                            ly[i] + ml_k_a,
                            lz[i] + ml_m_a,
                            lx[j] + r_x_c - ml_i_b,
                            ly[j] + r_y_c - ml_k_b,
                            lz[j] + r_z_c - ml_m_b,
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
                            lx[i] + ml_i_a,
                            ly[i] + ml_k_a,
                            lz[i] + ml_m_a,
                            lx[j] - ml_i_b,
                            ly[j] - ml_k_b,
                            lz[j] - ml_m_b,
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
                        - 2.0*exp[j]*ml_a[i]*(
                            nuclear_attraction(
                            lx[i] - ml_i_a,
                            ly[i] - ml_k_a,
                            lz[i] - ml_m_a,
                            lx[j] + r_x_c + ml_i_b,
                            ly[j] + r_y_c + ml_k_b,
                            lz[j] + r_z_c + ml_m_b,
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
                            lx[i] - ml_i_a,
                            ly[i] - ml_k_a,
                            lz[i] - ml_m_a,
                            lx[j] + ml_i_b,
                            ly[j] + ml_k_b,
                            lz[j] + ml_m_b,
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
                        + ml_a[i]*ml_b[j]*(
                            nuclear_attraction(
                            lx[i] - ml_i_a,
                            ly[i] - ml_k_a,
                            lz[i] - ml_m_a,
                            lx[j] + r_x_c - ml_i_b,
                            ly[j] + r_y_c - ml_k_b,
                            lz[j] + r_z_c - ml_m_b,
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
                            lx[i] - ml_i_a,
                            ly[i] - ml_k_a,
                            lz[i] - ml_m_a,
                            lx[j] - ml_i_b,
                            ly[j] - ml_k_b,
                            lz[j] - ml_m_b,
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
                )

            pnstcgop[count] = (
                -1.0 # p * A^2 * p = i*nabla*phi A^2 *i*nabla*phi => (i*i = -1)
                * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
                * (nefab + nefc)
                * 0.5
            )
            count += 1

    if output > 10:
        driver_time.add_name_delta_time(
            name = f"Gradient od diamagnetic Nuclear Shielding Tensor Atomic Integrals \
                for {magnetic_component} Magnetic Component and {spatial_sym} Spatial Symmetry",
                delta_time = (time() - start)
        )

    return pnstcgop

if __name__ == "__main__":
    # 6-311++G**
    pa2p = pnstcgop(
        coord = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        gauge = [0.0, 0.0, 0.0],
        spatial_sym=0,
        magnetic_component=0,
        atom=0,
        exp = [
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
        center = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        lx = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        ly = [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        lz = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        output=1,
        dalton_normalization = False,
        driver_time = None)

    print("P NSTGO P : ",pa2p)