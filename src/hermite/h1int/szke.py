from lib1h import *
from second_derivatives import *

############# Calculate the potential one body integrals ########################
def szke(
    coord, component, exp, center, lx, ly, lz, output, dalton_normalization, driver_time
):
    """_summary_

    Calculate the operator relationship with spin zeeman kinetic 
    derivative in the LRESC framework

            Nabla . Nabla^T + Nabla^T . Nabla

    Args:
        coord (list): list 2d with coordinates of the atoms
        component (int): Component to take of double derivative
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
        output (int): Output level for integral calculation
        dalton_normalization (bool): it is used the dalton normalization formule
        drive_time (drv_object): Object to manage the time

    Return:
        szke (array): array 1d with atomic integrals
    """

    start = time()
    if component == 0: # 2dxdx + dydy + dzdz
        dxdx: list = second_derivatives(coord, 0, exp, center, 
                            lx, ly, lz, output, dalton_normalization, driver_time)
        dydy: list = second_derivatives(coord, 1, exp, center, 
                                   lx, ly, lz, output, dalton_normalization, driver_time)
        dzdz: list = second_derivatives(coord, 2, exp, center, 
                                   lx, ly, lz, output, dalton_normalization, driver_time)
        szke: list = [2.0*a + b + c for a, b, c in zip(dxdx, dydy, dzdz)]
    elif component == 1 or component == 3: # dxdy or dydx
        szke = second_derivatives(coord, 3, exp, center, 
                                  lx, ly, lz, output, dalton_normalization, driver_time)
    elif component == 2 or component == 6: # dxdz or dzdx
        szke = second_derivatives(coord, 4, exp, center, 
                                  lx, ly, lz, output, dalton_normalization, driver_time)
    elif component == 4: # dxdx + 2dydy + dzdz
        dxdx = second_derivatives(coord, 0, exp, center, 
                            lx, ly, lz, output, dalton_normalization, driver_time)
        dydy = second_derivatives(coord, 1, exp, center, 
                                   lx, ly, lz, output, dalton_normalization, driver_time)
        dzdz = second_derivatives(coord, 2, exp, center, 
                                   lx, ly, lz, output, dalton_normalization, driver_time)
        szke = [a + 2.0*b + c for a, b, c in zip(dxdx, dydy, dzdz)]
    elif component == 5 or component == 7: # dydz or dzdy
        szke = second_derivatives(coord, 6, exp, center, 
                                  lx, ly, lz, output, dalton_normalization, driver_time)
    elif component == 8: # dxdx + dydy + 2dzdz
        dxdx = second_derivatives(coord, 0, exp, center, 
                            lx, ly, lz, output, dalton_normalization, driver_time)
        dydy = second_derivatives(coord, 1, exp, center, 
                                   lx, ly, lz, output, dalton_normalization, driver_time)
        dzdz = second_derivatives(coord, 2, exp, center, 
                                   lx, ly, lz, output, dalton_normalization, driver_time)
        szke = [a + b + 2.0*c for a, b, c in zip(dxdx, dydy, dzdz)]
    else:
        raise  ValueError(f"***Error\n\n This term of Nabla . Nabla^T + Nabla^T . Nabla not exist: {component}")

    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Laplacian component {component}", delta_time=(time() - start)
        )

    return szke

if __name__ == "__main__":
    # STO-2G
    print("\n LiH \n")
    s = szke(
        coord=[[0.0, 0.0, -0.545857052],[0.0, 0.0, 2.309057052]],
        component=2,
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
    print("szke : ", s, "\n", len(s), "\n\n")