from lib1h import *


def ppsop(
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

    ppsop: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    count: int = 0

    """
    Component Selection L = p x r
                        = (ypz-zpy)x + (zpx-xpz)y + (xpy-ypx)z
    where r = r_e - r_k
    """
    # Constants relationed with PSO operator
    r_x_left: int = 0
    r_y_left: int = 0
    r_z_left: int = 0
    r_x_right: int = 0
    r_y_right: int = 0
    r_z_right: int = 0

    if spatial_sym == 0:
        """X Component"""
        r_y_left = 1
        r_z_right =  1

        ml_i: int = lx
        ml_j: int = ly
        ml_l: int = lz

        pi: list = [(-1,0,0,-1,0,-1),( 1,0,0,-1,0,-1),
                    ( 1,0,0,-1,0, 1),(-1,0,0,-1,0, 1),
                    (-1,0,0, 1,0,-1),( 1,0,0, 1,0,-1),
                    (-1,0,0, 1,0, 1),( 1,0,0, 1,0, 1)]
        pj: list = [(0,-1,0,0,-1,-1),(0, 1,0,0,-1,-1),
                    (0, 1,0,0,-1, 1),(0,-1,0,0,-1, 1),
                    (0,-1,0,0, 1,-1),(0, 1,0,0, 1,-1),
                    (0,-1,0,0, 1, 1),(0, 1,0,0, 1, 1)]
        pl: list = [(0,0,-1,0,0,-2),(0,0,-1,0,0, 0),
                    (0,0,-1,0,0, 2),(0,0, 1,0,0,-2),
                    (0,0, 1,0,0, 0),(0,0, 1,0,0, 2)]
        pa: list = [(-1,0,0,-1,-1,0),( 1,0,0,-1,-1,0),
                    ( 1,0,0,-1, 1,0),(-1,0,0,-1, 1,0),
                    (-1,0,0, 1,-1,0),( 1,0,0, 1,-1,0),
                    (-1,0,0, 1, 1,0),( 1,0,0, 1, 1,0)]
        pb: list = [(0,0,-1,0,-1,-1),(0,0, 1,0,-1,-1),
                    (0,0, 1,0, 1,-1),(0,0,-1,0, 1,-1),
                    (0,0,-1,0,-1, 1),(0,0, 1,0,-1, 1),
                    (0,0,-1,0, 1, 1),(0,0, 1,0, 1, 1)]
        pc: list = [(0,-1,0,0,-2,0),(0,-1,0,0, 0,0),
                    (0,-1,0,0, 2,0),(0, 1,0,0,-2,0),
                    (0, 1,0,0, 0,0),(0, 1,0,0, 2,0)]
    elif spatial_sym == 1:
        """Y Componente"""
        r_z_left = 1
        r_x_right = 1

        # zpx - xpz
        ml_i: int = ly
        ml_j: int = lz
        ml_l: int = lx
        pi: list = [(0,-1,0,-1,-1,0),(0,1,0,-1,-1,0),
                    (0,1,0,1,-1,0),(0,-1,0,1,-1,0),
                    (0,-1,0,-1,1,0),(0,1,0,-1,1,0),
                    (0,-1,0,1,1,0),(0,1,0,1,1,0)]
        pj: list = [(0,0,-1,-1,0,-1),(0,0,1,-1,0,-1),
                    (0,0,1,1,0,-1),(0,0,-1,1,0,-1),
                    (0,0,-1,-1,0,1),(0,0,1,-1,0,1),
                    (0,0,-1,1,0,1),(0,0,1,1,0,1)]
        pl: list = [(-1,0,0,-2,0,0),(-1,0,0,0,0,0),
                    (-1,0,0,2,0,0),(1,0,0,-2,0,0),
                    (1,0,0,0,0,0),(1,0,0,2,0,0)]
        pa: list = [(0,-1,0,0,-1,-1),(0, 1,0,0,-1,-1),
                    (0, 1,0,0,-1, 1),(0,-1,0,0,-1, 1),
                    (0,-1,0,0, 1,-1),(0, 1,0,0, 1,-1),
                    (0,-1,0,0, 1, 1),(0, 1,0,0, 1, 1)]
        pb: list = [(-1,0,0,-1,0,-1),( 1,0,0,-1,0,-1),
                    ( 1,0,0,-1,0, 1),(-1,0,0,-1,0, 1),
                    (-1,0,0, 1,0,-1),( 1,0,0, 1,0,-1),
                    (-1,0,0, 1,0, 1),( 1,0,0, 1,0, 1)]
        pc: list = [(0,0,-1,0,0,-2),(0,0,-1,0,0, 0),
                    (0,0,-1,0,0, 2),(0,0, 1,0,0,-2),
                    (0,0, 1,0,0, 0),(0,0, 1,0,0, 2)]


    elif spatial_sym == 2:
        """Z Componente"""
        r_x_left = 1
        r_y_right = 1

        # xpy - ypx
        ml_i: int = lz
        ml_j: int = lx
        ml_l: int = ly

        pi: list = [(0,0,-1,0,-1,-1),(0,0,1,0,-1,-1),
                    (0,0,1,0,1,-1),(0,0,-1,0,1,-1),
                    (0,0,-1,0,-1,1),(0,0,1,0,-1,1),
                    (0,0,-1,0,1,1),(0,0,1,0,1,1)]
        pj: list = [(-1,0,0,-1,-1,0),(1,0,0,-1,-1,0),
                    (1,0,0,-1,1,0),(-1,0,0,-1,1,0),
                    (-1,0,0,1,-1,0),(1,0,0,1,-1,0),
                    (-1,0,0,1,1,0),(1,0,0,1,1,0)]
        pl: list = [(0,-1,0,0,-2,0),(0,-1,0,0,0,0),
                    (0,-1,0,0,2,0),(0,1,0,0,-2,0),
                    (0,1,0,0,0,0),(0,1,0,0,2,0)]
        pa: list = [(0,0,-1,-1,0,-1),(0,0,1,-1,0,-1),
                    (0,0,1,1,0,-1),(0,0,-1,1,0,-1),
                    (0,0,-1,-1,0,1),(0,0,1,-1,0,1),
                    (0,0,-1,1,0,1),(0,0,1,1,0,1)]
        pb: list = [(0,-1,0,-1,-1,0),(0,1,0,-1,-1,0),
                    (0,1,0,1,-1,0),(0,-1,0,1,-1,0),
                    (0,-1,0,-1,1,0),(0,1,0,-1,1,0),
                    (0,-1,0,1,1,0),(0,1,0,1,1,0)]
        pc: list = [(-1,0,0,-2,0,0),(-1,0,0,0,0,0),
                    (-1,0,0,2,0,0),(1,0,0,-2,0,0),
                    (1,0,0,0,0,0),(1,0,0,2,0,0)]
    else:
        raise ValueError(f"***Error\n\n Component not exist: {spatial_sym}")


    for i in range(total_nprim):

        for j in range(i, total_nprim):

            # left 
            ctes_left_ipsomi: list = [
                -ml_i[i]*ml_i[j]*ml_l[j],2.*exp[i]*ml_i[j]*ml_l[j],
                -4.*exp[i]*exp[j]*ml_i[j],2.*exp[j]*ml_i[i]*ml_i[j],
                2.*exp[j]*ml_i[i]*ml_l[j],-4.*exp[i]*exp[j]*ml_l[j],
                -4.*exp[j]*exp[j]*ml_i[i],8.*exp[i]*exp[j]*exp[j]
            ]
            ctes_left_jpsomj: list = [
                -ml_j[i]*ml_j[j]*ml_l[j],2.*exp[i]*ml_j[j]*ml_l[j],
                -4.*exp[i]*exp[j]*ml_j[j],2.*exp[j]*ml_j[i]*ml_j[j],
                2.*exp[j]*ml_j[i]*ml_l[j],-4.*exp[i]*exp[j]*ml_l[j],
                -4.*exp[j]*exp[j]*ml_j[i],8.*exp[i]*exp[j]*exp[j]
            ]
            ctes_left_lpsoml = [
                -ml_l[i]*(ml_l[j]*(ml_l[j]-1)),2.0*ml_l[i]*exp[j]*(2.0*ml_l[j]+1.0),
                -4.0*ml_l[i]*exp[j]*exp[j],2.0*exp[i]*(ml_l[j]*(ml_l[j]-1)),
                -4.0*exp[i]*exp[j]*(2.0*ml_l[j]+1.0),8.0*exp[i]*exp[j]*exp[j]
            ]
            # right
            ctes_right_ipsomi: list = [
                -ml_i[i]*ml_i[j]*ml_j[j],2.*exp[i]*ml_i[j]*ml_j[j],
                -4.*exp[i]*exp[j]*ml_i[j],2.*exp[j]*ml_i[i]*ml_i[j],
                2.*exp[j]*ml_i[i]*ml_j[j],-4.*exp[i]*exp[j]*ml_j[j],
                -4.*exp[j]*exp[j]*ml_i[i],8.*exp[i]*exp[j]*exp[j]
            ]
            ctes_right_jpsomj: list = [
                -ml_l[i]*ml_l[j]*ml_j[j],2.*exp[i]*ml_l[j]*ml_j[j],
                -4.*exp[i]*exp[j]*ml_l[j],2.*exp[j]*ml_l[i]*ml_l[j],
                2.*exp[j]*ml_l[i]*ml_j[j],-4.*exp[i]*exp[j]*ml_j[j],
                -4.*exp[j]*exp[j]*ml_l[i],8.*exp[i]*exp[j]*exp[j]
            ]
            ctes_right_lpsoml = [
                -ml_j[i]*(ml_j[j]*(ml_j[j]-1)),2.0*ml_j[i]*exp[j]*(2.0*ml_j[j]+1.0),
                -4.0*ml_j[i]*exp[j]*exp[j],2.0*exp[i]*(ml_j[j]*(ml_j[j]-1)),
                -4.0*exp[i]*exp[j]*(2.0*ml_j[j]+1.0),8.0*exp[i]*exp[j]*exp[j]
            ]
            countd: int = 0
            i_left_psoj_i: float = 0.0
            j_left_psoj_j: float = 0.0
            l_left_psoj_l: float = 0.0
            i_right_psoj_i: float = 0.0
            j_right_psoj_j: float = 0.0
            l_right_psoj_l: float = 0.0
            for ctes in ctes_left_ipsomi:
                # LEFT
                i_left_psoj_i -= ctes*nuclear_attraction(
                lx[i] + pi[countd][0],
                ly[i] + pi[countd][1],
                lz[i] + pi[countd][2],
                lx[j] + pi[countd][3],
                ly[j] + pi[countd][4],
                lz[j] + pi[countd][5],
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
                j_left_psoj_j -= ctes_left_jpsomj[countd]*nuclear_attraction(
                lx[i] + pj[countd][0],
                ly[i] + pj[countd][1],
                lz[i] + pj[countd][2],
                lx[j] + pj[countd][3],
                ly[j] + pj[countd][4],
                lz[j] + pj[countd][5],
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
                if countd < 6:
                    l_left_psoj_l -= ctes_left_lpsoml[countd]*nuclear_attraction(
                    lx[i] + pl[countd][0],
                    ly[i] + pl[countd][1],
                    lz[i] + pl[countd][2],
                    lx[j] + pl[countd][3],
                    ly[j] + pl[countd][4],
                    lz[j] + pl[countd][5],
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
            
                # RIGHT
                i_right_psoj_i -= ctes_right_ipsomi[countd]*nuclear_attraction(
                lx[i] + pa[countd][0],
                ly[i] + pa[countd][1],
                lz[i] + pa[countd][2],
                lx[j] + pa[countd][3],
                ly[j] + pa[countd][4],
                lz[j] + pa[countd][5],
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
                j_right_psoj_j -= ctes_right_jpsomj[countd]*nuclear_attraction(
                lx[i] + pb[countd][0],
                ly[i] + pb[countd][1],
                lz[i] + pb[countd][2],
                lx[j] + pb[countd][3],
                ly[j] + pb[countd][4],
                lz[j] + pb[countd][5],
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
                if countd < 6:
                    l_right_psoj_l -= ctes_right_lpsoml[countd]*nuclear_attraction(
                    lx[i] + pc[countd][0],
                    ly[i] + pc[countd][1],
                    lz[i] + pc[countd][2],
                    lx[j] + pc[countd][3],
                    ly[j] + pc[countd][4],
                    lz[j] + pc[countd][5],
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
                countd += 1



            ppsop[count] = (
                -normalization(lx[i], ly[i], lz[i], exp[i], dalton_normalization)
                * normalization(lx[j], ly[j], lz[j], exp[j], dalton_normalization)
                * 2.0 * np.pi / (exp[i] + exp[j])
                * (
                    (i_left_psoj_i  - i_right_psoj_i)
                    +(j_left_psoj_j - l_right_psoj_l)
                    +(l_left_psoj_l - j_right_psoj_j)
                )
            )
            count += 1


    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Paramagnetic Spin-Orbit Atomic Integrals, \
        for {spatial_sym} Spatial Symmetry",
            delta_time=(time() - start),
        )

    return ppsop

if __name__ == "__main__":
    # 6-311++G**
    # s = ppsop(
    #     coord=[[0.0, 0.0, -0.545857052],[0.0, 0.0, 2.309057052]],
    #     exp=[
    #         33.865,
    #         5.09479,
    #         1.15879,
    #         0.32584,
    #         0.102741,
    #         0.036,
    #         0.75,
    #         0.75,
    #         0.75,
    #         33.865,
    #         5.09479,
    #         1.15879,
    #         0.32584,
    #         0.102741,
    #         0.036,
    #         0.75,
    #         0.75,
    #         0.75,
    #     ],
    #     center=[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    #     lx=[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    #     ly=[0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    #     lz=[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    #     spatial_sym=0,
    #     atom=0,
    #     output=9,
    #     dalton_normalization=False,
    #     driver_time=None,
    # )
    
    # STO-2G
    s = ppsop(
        coord=[[0.0, 0.0, -0.545857052],[0.0, 0.0, 2.309057052]],
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
        spatial_sym=0,
        atom=0,
        output=9,
        dalton_normalization=False,
        driver_time=None,
    )

    print("ppsop : ", s, len(s))