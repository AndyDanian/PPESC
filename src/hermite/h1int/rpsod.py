from lib1h import *


def rpsod(
    coord,
    spatial_sym,
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
    Angular moment integrals by dipole and derivaties

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

    rpsod: list = [0 for i in range(int(total_nprim * total_nprim))]

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

    d_x: int = 0
    d_y: int = 0
    d_z: int = 0
    r_x: int = 0
    r_y: int = 0
    r_z: int = 0

    d2_x_left: int = 0
    d2_y_left: int = 0
    d2_z_left: int = 0
    d2_x_right: int = 0
    d2_y_right: int = 0
    d2_z_right: int = 0

    if spatial_sym == 0:
        """X Component"""
        # d_x = 1
        der_y_right = 1
        der_z_left = 1
        r_x = 1
        r_y_left = 1
        r_z_right = 1
        mli: list = lx
        l_right: list = ly
        l_left: list = lz
        d2_y_left = 2
        d2_z_right = 2
    elif spatial_sym == 1:
        """Y Componente"""
        d_y = 1
        der_z_right = 1
        der_x_left = 1
        r_y = 1
        r_x_right = 1
        r_z_left = 1
        mli: list = ly
        l_right: list = lz
        l_left: list = lx
        d2_z_left = 2
        d2_x_right = 2
    elif spatial_sym == 2:
        """Z Componente"""
        d_z = 1
        der_x_right = 1
        der_y_left = 1
        r_z = 1
        r_y_right = 1
        r_x_left = 1
        mli: list = lz
        l_right: list = lx
        l_left: list = ly
        d2_x_left = 2
        d2_y_right = 2
    else:
        raise ValueError(f"***Error\n\n Component not exist: {spatial_sym}")

    for i in range(total_nprim):

        for j in range(total_nprim):

            left_term: float = 0.0
            right_term: float = 0.0
            for d_x, d_y, d_z in zip([1, 0, 0],[0, 1, 0],[0, 0, 1]):
                if (der_x_left + d_x > 1 or der_y_left + d_y > 1 or der_z_left + d_z > 1):
                    left_term += (4.0 * exp[j] *  exp[j] *
                                nuclear_attraction(
                                lx[i],
                                ly[i],
                                lz[i],
                                lx[j] + der_x_left + d_x,
                                ly[j] + der_y_left + d_y,
                                lz[j] + der_z_left + d_z,
                                r_x_left + d_x,
                                r_y_left + d_y,
                                r_z_left + d_z,
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
                                - 2.0 * exp[j] * (2.0 * l_left[j] + 1.0) *
                                nuclear_attraction(
                                lx[i],
                                ly[i],
                                lz[i],
                                lx[j],
                                ly[j],
                                lz[j],
                                r_x_left + d_x,
                                r_y_left + d_y,
                                r_z_left + d_z,
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
                            + l_left[j] * (l_left[j] - 1.0) *
                                nuclear_attraction(
                                lx[i],
                                ly[i],
                                lz[i],
                                lx[j] - der_x_left - d_x,
                                ly[j] - der_y_left - d_y,
                                lz[j] - der_z_left - d_z,
                                r_x_left + d_x,
                                r_y_left + d_y,
                                r_z_left + d_z,
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

                    right_term += (
                                2.0 * exp[j] *
                                (
                                2.0 * exp[j] *
                                nuclear_attraction(
                                lx[i],
                                ly[i],
                                lz[i],
                                lx[j] + der_x_right + d_x,
                                ly[j] + der_y_right + d_y,
                                lz[j] + der_z_right + d_z,
                                r_x_right + d_x,
                                r_y_right + d_y,
                                r_z_right + d_z,
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
                            -
                            l_right[j] * 
                                nuclear_attraction(
                                lx[i],
                                ly[i],
                                lz[i],
                                lx[j] - der_x_right + d_x,
                                ly[j] - der_y_right + d_y,
                                lz[j] - der_z_right + d_z,
                                r_x_right + d_x,
                                r_y_right + d_y,
                                r_z_right + d_z,
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
                                - l_left[j] * 
                                (
                                    2.0 * exp[j] *
                                        nuclear_attraction(
                                        lx[i],
                                        ly[i],
                                        lz[i],
                                        lx[j] + der_x_right - d_x,
                                        ly[j] + der_y_right - d_y,
                                        lz[j] + der_z_right - d_z,
                                        r_x_right + d_x,
                                        r_y_right + d_y,
                                        r_z_right + d_z,
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
                                    ) - l_right[j] * 
                                    nuclear_attraction(
                                        lx[i],
                                        ly[i],
                                        lz[i],
                                        lx[j] - der_x_right - d_x,
                                        ly[j] - der_y_right - d_y,
                                        lz[j] - der_z_right - d_z,
                                        r_x_right + d_x,
                                        r_y_right + d_y,
                                        r_z_right + d_z,
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
                                )
                    # if count == 140:
                    #     print("Primer IF")
                    #     print("der_i_left ",der_x_left,der_y_left,der_z_left)
                    #     print("r_i_left ",r_x_left,r_y_left,r_z_left)
                    #     print("der_i_right ",der_x_right,der_y_right,der_z_right)
                    #     print("r_i_right ",r_x_right,r_y_right,r_z_right)
                    #     print("d_i ",d_x,d_y,d_z)
                    #     print(" LEFT TERM : ",left_term * 2.0
                    #                         * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                    #                         * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                    #                         * 2.0
                    #                         * np.pi
                    #                         / (exp[i] + exp[j]))
                    #     print("l_right ",l_right[j]," lz ",lz[j])
                    #     print(" RIGHT TERM : ",right_term * 2.0
                    #                         * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                    #                         * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                    #                         * 2.0
                    #                         * np.pi
                    #                         / (exp[i] + exp[j]))
                    #     print(right_term)
                    #     print()
                elif (der_x_right + d_x > 1 or der_y_right + d_y > 1 or der_z_right + d_z > 1):
                    right_term += (4.0 * exp[j] *  exp[j] *
                                nuclear_attraction(
                                lx[i],
                                ly[i],
                                lz[i],
                                lx[j] + der_x_right + d_x,
                                ly[j] + der_y_right + d_y,
                                lz[j] + der_z_right + d_z,
                                r_x_right + d_x,
                                r_y_right + d_y,
                                r_z_right + d_z,
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
                                - 2.0 * exp[j] * (2.0 * l_right[j] + 1.0) *
                                nuclear_attraction(
                                lx[i],
                                ly[i],
                                lz[i],
                                lx[j],
                                ly[j],
                                lz[j],
                                r_x_right + d_x,
                                r_y_right + d_y,
                                r_z_right + d_z,
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
                            + l_right[j] * (l_right[j] - 1.0) *
                                nuclear_attraction(
                                lx[i],
                                ly[i],
                                lz[i],
                                lx[j] - der_x_right - d_x,
                                ly[j] - der_y_right - d_y,
                                lz[j] - der_z_right - d_z,
                                r_x_right + d_x,
                                r_y_right + d_y,
                                r_z_right + d_z,
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

                    left_term += (2.0 * exp[j] * (
                        2.0 * exp[j] * 
                            nuclear_attraction(
                            lx[i],
                            ly[i],
                            lz[i],
                            lx[j] + der_x_left + d_x,
                            ly[j] + der_y_left + d_y,
                            lz[j] + der_z_left + d_z,
                            r_x_left + d_x,
                            r_y_left + d_y,
                            r_z_left + d_z,
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
                            lx[j] - der_x_left + d_x,
                            ly[j] - der_y_left + d_y,
                            lz[j] - der_z_left + d_z,
                            r_x_left + d_x,
                            r_y_left + d_y,
                            r_z_left + d_z,
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
                        ))-l_right[j]*(
                            2.0 * exp[j] * nuclear_attraction(
                                lx[i],
                                ly[i],
                                lz[i],
                                lx[j] + der_x_left - d_x,
                                ly[j] + der_y_left - d_y,
                                lz[j] + der_z_left - d_z,
                                r_x_left + r_x,
                                r_y_left + r_y,
                                r_z_left + r_z,
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
                                lx[j] - der_x_left - d_x,
                                ly[j] - der_y_left - d_y,
                                lz[j] - der_z_left - d_z,
                                r_x_left + r_x,
                                r_y_left + r_y,
                                r_z_left + r_z,
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
                        )
                    # if count == 140:
                    #     print("Segundo IF")
                    #     print("der_i_right ",der_x_right,der_y_right,der_z_right)
                    #     print("r_i_right ",r_x_right,r_y_right,r_z_right)
                    #     print("d_i ",d_x,d_y,d_z)
                    #     print(" LEFT TERM : ",left_term * 2.0
                    #                         * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                    #                         * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                    #                         * 2.0
                    #                         * np.pi
                    #                         / (exp[i] + exp[j]))
                    #     print(" RIGHT TERM : ",right_term * 2.0
                    #                         * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                    #                         * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                    #                         * 2.0
                    #                         * np.pi
                    #                         / (exp[i] + exp[j]))
                    #     print()
                else:
                    left_term += (2.0 * exp[j] * (
                        2.0 * exp[j] * 
                            nuclear_attraction(
                            lx[i],
                            ly[i],
                            lz[i],
                            lx[j] + der_x_left + d_x,
                            ly[j] + der_y_left + d_y,
                            lz[j] + der_z_left + d_z,
                            r_x_left + d_x,
                            r_y_left + d_y,
                            r_z_left + d_z,
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
                            lx[j] - der_x_left + d_x,
                            ly[j] - der_y_left + d_y,
                            lz[j] - der_z_left + d_z,
                            r_x_left + d_x,
                            r_y_left + d_y,
                            r_z_left + d_z,
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
                        ))-mli[j]*(
                            2.0 * exp[j] * nuclear_attraction(
                                lx[i],
                                ly[i],
                                lz[i],
                                lx[j] + der_x_left - d_x,
                                ly[j] + der_y_left - d_y,
                                lz[j] + der_z_left - d_z,
                                r_x_left + r_x,
                                r_y_left + r_y,
                                r_z_left + r_z,
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
                                lx[j] - der_x_left - d_x,
                                ly[j] - der_y_left - d_y,
                                lz[j] - der_z_left - d_z,
                                r_x_left + r_x,
                                r_y_left + r_y,
                                r_z_left + r_z,
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
                        )

                    right_term += (2.0 * exp[j] * (
                                2.0 * exp[j] *
                                nuclear_attraction(
                                lx[i],
                                ly[i],
                                lz[i],
                                lx[j] + der_x_right + d_x,
                                ly[j] + der_y_right + d_y,
                                lz[j] + der_z_right + d_z,
                                r_x_right + d_x,
                                r_y_right + d_y,
                                r_z_right + d_z,
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
                                lx[j] - der_x_right + d_x,
                                ly[j] - der_y_right + d_y,
                                lz[j] - der_z_right + d_z,
                                r_x_right + d_x,
                                r_y_right + d_y,
                                r_z_right + d_z,
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
                                )) - mli[j] * (
                                    2.0 * exp[j] *
                                        nuclear_attraction(
                                        lx[i],
                                        ly[i],
                                        lz[i],
                                        lx[j] + der_x_right - d_x,
                                        ly[j] + der_y_right - d_y,
                                        lz[j] + der_z_right - d_z,
                                        r_x_right + d_x,
                                        r_y_right + d_y,
                                        r_z_right + d_z,
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
                                        lx[j] - der_x_right - d_x,
                                        ly[j] - der_y_right - d_y,
                                        lz[j] - der_z_right - d_z,
                                        r_x_right + d_x,
                                        r_y_right + d_y,
                                        r_z_right + d_z,
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
                                )

                    # if count == 140:
                    #     print("ELSE")
                    #     print()
                    #     print("der_i_left ",der_x_left,der_y_left,der_z_left)
                    #     print("r_i_left ",r_x_left,r_y_left,r_z_left)
                    #     print("d_i ",d_x,d_y,d_z)
                    #     print()
                    #     print("der_i_right ",der_x_right,der_y_right,der_z_right)
                    #     print("r_i_right ",r_x_right,r_y_right,r_z_right)
                    #     print("d_i ",d_x,d_y,d_z)
                    #     print(" LEFT TERM : ",left_term * 2.0
                    #                         * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                    #                         * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                    #                         * 2.0
                    #                         * np.pi
                    #                         / (exp[i] + exp[j]))
                    #     print(" RIGHT TERM : ",right_term * 2.0
                    #                         * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                    #                         * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                    #                         * 2.0
                    #                         * np.pi
                    #                         / (exp[i] + exp[j]))
                    #     print()

            rpsod[count] = (
                2.0
                * normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
                * (left_term - right_term)
            )
            #!? WARNING los terminos z_k²dydz/r⁵_k dan con una diferencia a Wolfram Alpha.
            #!? en Wolfram estos salen con Warnings 
            # if count == 140:
            if abs(rpsod[count]) > 0.0:
                print("alpha ",exp[i], lx[i], ly[i], lz[i])
                print("beta  ",exp[j], lx[j], ly[j], lz[j])
                print(count," : ",rpsod[count])
                print()
            count += 1
            
    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Paramagnetic Spin-Orbit Atomic Integrals By Dipole and Derivaties, \
        for {spatial_sym} Spatial Symmetry",
            delta_time=(time() - start),
        )

    return rpsod

if __name__ == "__main__":
    # STO-2G
    print("\n LiH \n")
    s = rpsod(
        spatial_sym=0,
        coord=[[0.0, 0.0, -0.545857052],[0.0, 0.0, 2.309057052]],
        atom=0,
        exp=[
            6.1638450, 1.0971610, 0.2459160, 0.0623710,
            0.2459160, 0.2459160, 0.2459160,
            0.0623710, 0.0623710, 0.0623710,
            1.3097564, 0.2331360
        ],
        center=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
        lx=[0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0],
        ly=[0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0],
        lz=[0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0],
        output=9,
        dalton_normalization=False,
        driver_time=None,
    )
    print("rpsod : ", s, "\n", len(s), "\n\n")
