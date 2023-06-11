from lib1h import *


def psopso(
    coord,
    spatial_sym,
    magnetic_component,
    atom,
    exp,
    center,
    lx,
    ly,
    lz,
    output,
    dalton_normalization,
    driver_time,
):
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
        dalton_normalization (bool): it is used the dalton normalization formule
        drive_time (drv_object): Object to manage the time

    Return:
        psooz (array): array 1d with atomic integrals
    """

    start: float = time()
    # Primitive total in the cluster
    total_nprim: int = len(exp)

    psopso: list = [0 for i in range(int(total_nprim * total_nprim))]

    count: int = 0

    pso_a_r_x_b: int = 0
    pso_a_r_y_b: int = 0
    pso_a_r_z_b: int = 0
    pso_a_r_x_c: int = 0
    pso_a_r_y_c: int = 0
    pso_a_r_z_c: int = 0

    pso_b_r_x_b: int = 0
    pso_b_r_y_b: int = 0
    pso_b_r_z_b: int = 0
    pso_b_r_x_c: int = 0
    pso_b_r_y_c: int = 0
    pso_b_r_z_c: int = 0

    # * Associated with the Orbital Zemann Operator
    if magnetic_component == 0:
        """X Component"""
        pso_a_r_y_b = 1
        pso_a_r_z_c = 1
        pso_a_der_l_b: list = ly
        pso_a_der_l_c: list = lz
    elif magnetic_component == 1:
        """Y Component"""
        pso_a_r_z_b = 1
        pso_a_r_x_c = 1
        pso_a_der_l_b: list = lz
        pso_a_der_l_c: list = lx
    elif magnetic_component == 2:
        """Z Component"""
        pso_a_r_x_b = 1
        pso_a_r_y_c = 1
        pso_a_der_l_b: list = lx
        pso_a_der_l_c: list = ly

    # * Associated with the Paramagnetic Spin--Orbit Operator
    if spatial_sym == 0:
        """1 Component"""
        pso_b_r_y_b = 1
        pso_b_r_z_c = 1
        pso_b_der_l_b: list = ly
        pso_b_der_l_c: list = lz
    elif spatial_sym == 1:
        """2 Componente"""
        pso_b_r_z_b = 1
        pso_b_r_x_c = 1
        pso_b_der_l_b: list = lz
        pso_b_der_l_c: list = lx
    elif spatial_sym == 2:
        """3 Componente"""
        pso_b_r_x_b = 1
        pso_b_r_y_c = 1
        pso_b_der_l_b: list = lx
        pso_b_der_l_c: list = ly
    else:
        raise ValueError(f"***Error\n\n Component not exist: {spatial_sym}")

    for i in range(total_nprim):

        for j in range(total_nprim):
            """
            <phi_i|{L_k/r³_k, L_k/r³_k}|phi_j> : One operator is applied to the bra while other is applied a ket

            for one term of ati--commutator is

            f = x_a^iy_a^jz_a^kexp(-p r_a^2) (-i nabla x r_0) (-i nabla x r_k)/r³_k x_l^iy_m^jz_n^kexp(-t r_b^2)

            then

            f = -[
                {-r_a^i s_a^js_0 (2pt_a^{k+1}-kt_a^{k-1}) exp(-p r_a^2)
                +r_a^i (2ps_a^{j+1}-js_a^{j-1}) t_a^kt_0 exp(-p r_a^2)}
                {r_b^l s_b^ms_k/r³_k (2pt_n^{n+1}-nt_b^{n-1}) exp(-t r_b^2)
                -r_b^l (2ps_b^{m+1}-ms_b^{m-1}) t_b^nt_k/r³_k exp(-t r_b^2)}
                ]

            where r, s, t, u, v y w are x, y, or z. Also, the term in function of the gauge o k-th atom's coordinate
            are re-written

                                    x-0_x = (x-a/b_x) + (a/b_x-0_x)
            """

            igdj_jgdi_idj_jdi: float = 0.0
            igdj_jgdi_idj_jdi_b: float = 0.0
            for k in range(4):
                """
                                            (r_k d_s - s_k d_r)
                l_r_i_a  [int]: is the potentia or ml_i associated to the bra and position operator of the OZ (r)
                l_s_i_b  [int]: is the potentia or ml_i associated to the bra and derivator operator of the OZ (d)

                oz_r_t_b [int]: is the potentia of the position operator in OZ, associated with position operator (r)
                oz_s_t_c [int]: is the potentia of the position operator in OZ, associated with derivator operator (d)

                l_r_j    [int]: is the potentia or ml_i associated to the ket and derivator operator of the OZ (d)
                r_s      [int]: is the potentia associated with the operator d^i_s(1/r)

                                            d_iG = s^k_a r^l_a (2 p t_a^{m+1} - m t_a^{m-1}) exp(- p r²_a)
                pso_a_der_l_c  [int]: is the ml that result to apply d_i to the gaussian function (m) and associated with the bra
                                    and second derivative of OZ
                pso_a_der_l_b  [int]: is the ml that result to apply d_i to the gaussian function (m) and associated with the bra
                                    and first derivative of OZ
                pso_b_der_l_c [int]: is the ml that result to apply d_i to the gaussian function (m) and associated with the ket
                                    and second derivative of PSO
                pso_b_der_l_b [int]: is the ml that result to apply d_i to the gaussian function (m) and associated with the ket
                                    and first derivative of PSO


                atom        [int]: is indicated by the spatial symmetry or i.e. PSO
                """
                if k == 0:
                    """
                                            OZ.PSO
                    [r_a^i s_a^j (2pt_a^{k+1}-kt_a^{k-1}) s_k/r³_k exp(-p r_a^2)]
                    [r_b^l s_b^m (2up_b^{n+1}-np_b^{n-1}) s_k/r³_k exp(-u r_b^2)]
                                            PSO.OZ
                            <phi| r_k/r³_k d_s t_o d_u |phi>
                    """
                    # * PSO A
                    l_x_i_pso_a, l_y_i_pso_a, l_z_i_pso_a = (
                        pso_a_r_x_c,
                        pso_a_r_y_c,
                        pso_a_r_z_c,
                    )  # ? derivative
                    r_x_pso_a, r_y_pso_a, r_z_pso_a = (
                        pso_a_r_x_b,
                        pso_a_r_y_b,
                        pso_a_r_z_b,
                    )  # ? position
                    der_l_pso_a_bra: float = float(pso_a_der_l_c[i])
                    der_l_pso_a_ket: float = float(pso_a_der_l_c[j])
                    # * PSO B
                    l_x_j_pso_b, l_y_j_pso_b, l_z_j_pso_b = (
                        pso_b_r_x_c,
                        pso_b_r_y_c,
                        pso_b_r_z_c,
                    )  # ? derivative
                    r_x_pso_b, r_y_pso_b, r_z_pso_b = (
                        pso_b_r_x_b,
                        pso_b_r_y_b,
                        pso_b_r_z_b,
                    )  # ? position
                    der_l_pso_b_bra: float = float(pso_b_der_l_c[i])
                    der_l_pso_b_ket: float = float(pso_b_der_l_c[j])
                    # * r_{a,i}r_{b,j}/(r_a^3r_b^3)
                    coef: float = 1.0 / 3.0
                elif k == 1:
                    """
                    -[r_a^i (2ps_a^{j+1}-js_a^{j-1}) t_a^{k} t_k/r³_k exp(-p r_a^2)]
                    [r_b^l s_b^m (2ut_n^{n+1}-nt_b^{n-1})    s_k/r³_k exp(-u r_b^2)]
                    """
                    # * PSO A
                    l_x_i_pso_a, l_y_i_pso_a, l_z_i_pso_a = (
                        pso_a_r_x_b,
                        pso_a_r_y_b,
                        pso_a_r_z_b,
                    )
                    r_x_pso_a, r_y_pso_a, r_z_pso_a = (
                        pso_a_r_x_c,
                        pso_a_r_y_c,
                        pso_a_r_z_c,
                    )
                    der_l_pso_a_bra: float = float(pso_a_der_l_b[i])
                    der_l_pso_a_ket: float = float(pso_a_der_l_b[j])
                    # * PSO B
                    l_x_j_pso_b, l_y_j_pso_b, l_z_j_pso_b = (
                        pso_b_r_x_c,
                        pso_b_r_y_c,
                        pso_b_r_z_c,
                    )
                    r_x_pso_b, r_y_pso_b, r_z_pso_b = (
                        pso_b_r_x_b,
                        pso_b_r_y_b,
                        pso_b_r_z_b,
                    )
                    der_l_pso_b_bra: float = float(pso_b_der_l_c[i])
                    der_l_pso_b_ket: float = float(pso_b_der_l_c[j])
                    # * r_{a,i}r_{b,j}/(r_a^3r_b^3)
                    coef: float = -1.0 / 3.0
                elif k == 2:
                    """
                    -[r_a^i s_a^{j} (2pt_a^{k+1}-kt_a^{k-1}) s_k/r³_k exp(-p r_a^2)]
                    [r_b^l (2us_b^{m+1}-ms_b^{m-1}) t_b^n    t_k/r³_k exp(-u r_b^2)]
                    """
                    # * PSO A
                    l_x_i_pso_a, l_y_i_pso_a, l_z_i_pso_a = (
                        pso_a_r_x_c,
                        pso_a_r_y_c,
                        pso_a_r_z_c,
                    )
                    r_x_pso_a, r_y_pso_a, r_z_pso_a = (
                        pso_a_r_x_b,
                        pso_a_r_y_b,
                        pso_a_r_z_b,
                    )
                    der_l_pso_a_bra: float = float(pso_a_der_l_c[i])
                    der_l_pso_a_ket: float = float(pso_a_der_l_c[j])
                    # * PSO B
                    l_x_j_pso_b, l_y_j_pso_b, l_z_j_pso_b = (
                        pso_b_r_x_b,
                        pso_b_r_y_b,
                        pso_b_r_z_b,
                    )
                    r_x_pso_b, r_y_pso_b, r_z_pso_b = (
                        pso_b_r_x_c,
                        pso_b_r_y_c,
                        pso_b_r_z_c,
                    )
                    der_l_pso_b_bra: float = float(pso_b_der_l_b[i])
                    der_l_pso_b_ket: float = float(pso_b_der_l_b[j])
                    # * r_{a,i}r_{b,j}/(r_a^3r_b^3)
                    coef: float = -1.0 / 3.0
                elif k == 3:
                    """
                    [r_a^i (2ps_a^{j+1}-js_a^{j-1}) t_a^kt_k/r³_k exp(-p r_a^2)]
                    [r_b^l (2us_b^{m+1}-ms_b^{m-1}) t_b^nt_k/r³_k exp(-u r_b^2)]
                    """
                    # * PSO A
                    l_x_i_pso_a, l_y_i_pso_a, l_z_i_pso_a = (
                        pso_a_r_x_b,
                        pso_a_r_y_b,
                        pso_a_r_z_b,
                    )
                    r_x_pso_a, r_y_pso_a, r_z_pso_a = (
                        pso_a_r_x_c,
                        pso_a_r_y_c,
                        pso_a_r_z_c,
                    )
                    der_l_pso_a_bra: float = float(pso_a_der_l_b[i])
                    der_l_pso_a_ket: float = float(pso_a_der_l_b[j])
                    # * PSO B
                    l_x_j_pso_b, l_y_j_pso_b, l_z_j_pso_b = (
                        pso_b_r_x_b,
                        pso_b_r_y_b,
                        pso_b_r_z_b,
                    )
                    r_x_pso_b, r_y_pso_b, r_z_pso_b = (
                        pso_b_r_x_c,
                        pso_b_r_y_c,
                        pso_b_r_z_c,
                    )
                    der_l_pso_b_bra: float = float(pso_b_der_l_b[i])
                    der_l_pso_b_ket: float = float(pso_b_der_l_b[j])
                    # * r_{a,i}r_{b,j}/(r_a^3r_b^3)
                    coef: float = 1.0 / 3.0
                """
                                        OZ.PSO
                OZ is applied to the bra and PSO to the ket
                """
                igdj_jgdi_idj_jdi += coef * (
                    4.0
                    * exp[i]
                    * exp[j]
                    * nuclear_attraction(
                        # *  PSO A
                        lx[i] + l_x_i_pso_a,
                        ly[i] + l_y_i_pso_a,
                        lz[i] + l_z_i_pso_a,
                        # *
                        # * PSO B
                        lx[j] + l_x_j_pso_b,
                        ly[j] + l_y_j_pso_b,
                        lz[j] + l_z_j_pso_b,
                        # * PSO A + PSO B
                        r_x_pso_a + r_x_pso_b,
                        r_y_pso_a + r_y_pso_b,
                        r_z_pso_a + r_z_pso_b,
                        # *
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        # * PSO
                        coord[atom][0],
                        coord[atom][1],
                        coord[atom][2],
                        # *
                    )
                    - 2.0
                    * exp[i]
                    * der_l_pso_b_ket
                    * nuclear_attraction(
                        lx[i] + l_x_i_pso_a,
                        ly[i] + l_y_i_pso_a,
                        lz[i] + l_z_i_pso_a,
                        lx[j] - l_x_j_pso_b,
                        ly[j] - l_y_j_pso_b,
                        lz[j] - l_z_j_pso_b,
                        r_x_pso_a + r_x_pso_b,
                        r_y_pso_a + r_y_pso_b,
                        r_z_pso_a + r_z_pso_b,
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
                    * der_l_pso_a_bra
                    * nuclear_attraction(
                        lx[i] - l_x_i_pso_a,
                        ly[i] - l_y_i_pso_a,
                        lz[i] - l_z_i_pso_a,
                        lx[j] + l_x_j_pso_b,
                        ly[j] + l_y_j_pso_b,
                        lz[j] + l_z_j_pso_b,
                        r_x_pso_a + r_x_pso_b,
                        r_y_pso_a + r_y_pso_b,
                        r_z_pso_a + r_z_pso_b,
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
                    + der_l_pso_a_bra
                    * der_l_pso_b_ket
                    * nuclear_attraction(
                        lx[i] - l_x_i_pso_a,
                        ly[i] - l_y_i_pso_a,
                        lz[i] - l_z_i_pso_a,
                        lx[j] - l_x_j_pso_b,
                        ly[j] - l_y_j_pso_b,
                        lz[j] - l_z_j_pso_b,
                        r_x_pso_a + r_x_pso_b,
                        r_y_pso_a + r_y_pso_b,
                        r_z_pso_a + r_z_pso_b,
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
                """
                                        PSO.OZ
                PSO is applied to the bra while OZ to the ket
                """
                igdj_jgdi_idj_jdi_b += coef * (
                    4.0
                    * exp[i]
                    * exp[j]
                    * nuclear_attraction(
                        # *  PSO
                        lx[i] + l_x_j_pso_b,
                        ly[i] + l_y_j_pso_b,
                        lz[i] + l_z_j_pso_b,
                        # *
                        # * OZ
                        lx[j] + l_x_i_pso_a,
                        ly[j] + l_y_i_pso_a,
                        lz[j] + l_z_i_pso_a,
                        # * PSO
                        r_x_pso_a + r_x_pso_b,
                        r_y_pso_a + r_y_pso_b,
                        r_z_pso_a + r_z_pso_b,
                        # *
                        exp[i],
                        exp[j],
                        coord[center[i]][0],
                        coord[center[i]][1],
                        coord[center[i]][2],
                        coord[center[j]][0],
                        coord[center[j]][1],
                        coord[center[j]][2],
                        # * PSO
                        coord[atom][0],
                        coord[atom][1],
                        coord[atom][2],
                        # *
                    )
                    - 2.0
                    * exp[j]
                    * der_l_pso_b_bra
                    * nuclear_attraction(
                        lx[i] - l_x_j_pso_b,
                        ly[i] - l_y_j_pso_b,
                        lz[i] - l_z_j_pso_b,
                        lx[j] + l_x_i_pso_a,
                        ly[j] + l_y_i_pso_a,
                        lz[j] + l_z_i_pso_a,
                        r_x_pso_a + r_x_pso_b,
                        r_y_pso_a + r_y_pso_b,
                        r_z_pso_a + r_z_pso_b,
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
                    * der_l_pso_a_ket
                    * nuclear_attraction(
                        lx[i] + l_x_j_pso_b,
                        ly[i] + l_y_j_pso_b,
                        lz[i] + l_z_j_pso_b,
                        lx[j] - l_x_i_pso_a,
                        ly[j] - l_y_i_pso_a,
                        lz[j] - l_z_i_pso_a,
                        r_x_pso_a + r_x_pso_b,
                        r_y_pso_a + r_y_pso_b,
                        r_z_pso_a + r_z_pso_b,
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
                    + der_l_pso_a_ket
                    * der_l_pso_b_bra
                    * nuclear_attraction(
                        lx[i] - l_x_j_pso_b,
                        ly[i] - l_y_j_pso_b,
                        lz[i] - l_z_j_pso_b,
                        lx[j] - l_x_i_pso_a,
                        ly[j] - l_y_i_pso_a,
                        lz[j] - l_z_i_pso_a,
                        r_x_pso_a + r_x_pso_b,
                        r_y_pso_a + r_y_pso_b,
                        r_z_pso_a + r_z_pso_b,
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

                # if (i == 2 and j == 2) or (i == 0 and j ==0):
                #    print("\n ka ",k,i,j,igdj_jgdi_idj_jdi)

            psopso[count] = (
                normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
                * (igdj_jgdi_idj_jdi + igdj_jgdi_idj_jdi_b)
            )
            if abs(psopso[count]) > 0.01:
                print("  <phi|{PSO, OZ}|phi> : ", i + 1, j + 1, psopso[count])
            # normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
            # * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
            # *2.0
            # * np.pi
            # / (exp[i] + exp[j])
            # * igdj_jgdi_idj_jdi_b)
            # exit()
            count += 1

    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Orbital-Zeeman Correction to the Paramagnetic Spin-Orbit Atomic Integrals, \
        for {magnetic_component} Magnetic Component, {spatial_sym} Spatial Symmetry, and {atom}-th Atom",
            delta_time=(time() - start),
        )

    return psopso


if __name__ == "__main__":
    # STO-2G
    lih = False
    if lih:
        print("\n LiH \n")
        s = psopso(
            coord=[[0.0, 0.0, -0.545857052], [0.0, 0.0, 2.309057052]],
            exp=[
                6.1638450,
                1.0971610,
                0.2459160,
                0.0623710,
                0.2459160,
                0.2459160,
                0.2459160,
                0.0623710,
                0.0623710,
                0.0623710,
                1.3097564,
                0.2331360,
            ],
            center=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
            lx=[0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0],
            ly=[0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0],
            lz=[0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0],
            output=9,
            spatial_sym=0,
            magnetic_component=1,
            atom=0,
            dalton_normalization=False,
            driver_time=None,
        )
    else:
        print("\n He [1s1p1d] \n")
        s = psopso(
            coord=[[0, 0, 0]],
            exp=[
                9623.91395,
                6.25523565,
                6.25523565,
                6.25523565,
                4.32782104,
                4.32782104,
                4.32782104,
                4.32782104,
                4.32782104,
                4.32782104,
                2.68495795,
                2.68495795,
                2.68495795,
            ],
            center=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            lx=[0, 1, 0, 0, 2, 1, 1, 0, 0, 0, 3, 2, 2],
            ly=[0, 0, 1, 0, 0, 1, 0, 2, 1, 0, 0, 1, 0],
            lz=[0, 0, 0, 1, 0, 0, 1, 0, 1, 2, 0, 0, 1],
            output=9,
            spatial_sym=0,
            magnetic_component=0,
            atom=0,
            dalton_normalization=True,
            driver_time=None,
        )

    print("psopso : ", s, "\n", len(s), "\n\n")
