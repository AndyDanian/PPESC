from libh import *

def sd(coord, magnetic_component, spatial_sym, atom, exp, center, lx, ly, lz, output):
    """
    Spin dipolar atomic integrals, which is a vector

    Args:
        coord (list): list 2d with coordinates of the atoms
        magnetic_component (int): magnetic component
        spatial_sym (int): spatial symmetry index
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

    sd: list = [0 for i in range(int(total_nprim * (total_nprim + 1) / 2))]

    count: int = 0

    # ! Buscar ref gfactor
    GFACTOR: float = 2.0023193134
    CONST_SD: float = GFACTOR / 2.0 * 1 / 3.0

    dx_x: int = 0
    dy_x: int = 0
    dz_x: int = 0
    dx_y: int = 0
    dy_y: int = 0
    dz_y: int = 0
    dx_z: int = 0
    dy_z: int = 0
    dz_z: int = 0

    CTE_XR: float = 1.0
    CTE_YR: float = 1.0
    CTE_ZR: float = 1.0

    x_active: bool = False
    y_active: bool = False
    z_active: bool = False

    if magnetic_component == 0:
        if spatial_sym == 0:
            CTE_XR = 2.0
            CTE_YR = -1.0
            CTE_ZR = -1.0
            dx_x = 2
            dy_y = 2
            dz_z = 2
            x_active = True
            y_active = True
            z_active = True
        elif spatial_sym == 1:
            CTE_YR = 3.0
            y_active = True
            dx_y = 1
            dy_y = 1
        elif spatial_sym == 2:
            CTE_ZR = 3.0
            z_active = True
            dx_z = 1
            dz_z = 1
    if magnetic_component == 1:
        if spatial_sym == 0:
            CTE_XR = 3.0
            x_active = True
            dx_x = 1
            dy_x = 1
        elif spatial_sym == 1:
            CTE_XR = -1.0
            CTE_YR = 2.0
            CTE_ZR = -1.0
            dx_x = 2
            dy_y = 2
            dz_z = 2
            x_active = True
            y_active = True
            z_active = True
        elif spatial_sym == 2:
            CTE_ZR = 3.0
            z_active = True
            dy_z = 1
            dz_z = 1
    if magnetic_component == 2:
        if spatial_sym == 0:
            CTE_XR = 3.0
            x_active = True
            dx_x = 1
            dz_x = 1
        elif spatial_sym == 1:
            CTE_YR = 3.0
            y_active = True
            dy_y = 1
            dz_y = 1
        elif spatial_sym == 2:
            CTE_YR = -1.0
            CTE_XR = -1.0
            CTE_ZR = 2.0
            dx_x = 2
            dy_y = 2
            dz_z = 2
            x_active = True
            y_active = True
            z_active = True


    for i in range(total_nprim):

        for j in range(i, total_nprim):

            # SD = <phi|(3r_krk^T-r_k^2)/rk^5|phi>
            # https://www.wolframalpha.com/input?i2d=true&i=D%5BDivide%5B1%2CSqrt%5BPower%5Bx%2C2%5D%2BPower%5By%2C2%5D%2BPower%5Bz%2C2%5D%5D%5D%2C%7Bx%2C2%7D%5D
            # !Note: is add the constant 1/3, view the before wolfram result

            xr = yr = zr = 0.0
            if x_active:
                xr = nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j],
                ly[j],
                lz[j],
                dx_x,
                dy_x,
                dz_x,
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
            if y_active:
                yr = nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j],
                ly[j],
                lz[j],
                dx_y, 
                dy_y, 
                dz_y, 
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
            if z_active:
                zr = nuclear_attraction(
                lx[i],
                ly[i],
                lz[i],
                lx[j],
                ly[j],
                lz[j],
                dx_z,
                dy_z,
                dz_z,
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

            sd[count] = (
                CONST_SD
                * Norm[lx[i] + ly[i] + lz[i]](exp[i])
                * Norm[lx[j] + ly[j] + lz[j]](exp[j])
                * 2.0
                * np.pi
                / (exp[i] + exp[j])
                * (CTE_XR * xr + CTE_YR * yr + CTE_ZR * zr)
            )
            count += 1

    if output > 10:
        print(f"\n *** Spin dipolar atomic integrals,\
        spatial symmetry {spatial_sym} and magnetic component {magnetic_component}\
        of atom {atom + 1}-th, time [s]: {time() - start:.6f}")

    return sd
