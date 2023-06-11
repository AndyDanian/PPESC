from lib1h import *
from laplacian_didj import *

############# Calculate the potential one body integrals ########################
def szmv(
    coord, component, exp, center, lx, ly, lz, output, dalton_normalization, driver_time
):
    """_summary_

    Calculate the operator relationship with spin zeeman massvelocity
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
        szmv (array): array 1d with atomic integrals
    """

    start = time()
    if component == 0 or component == 4 or component == 8:  # 2dxdx + d2yd2y + d2zd2z
        d2xd2x: list = laplacian_didj(
            coord, 0, exp, center, lx, ly, lz, output, dalton_normalization, driver_time
        )
        d2xd2y: list = laplacian_didj(
            coord, 3, exp, center, lx, ly, lz, output, dalton_normalization, driver_time
        )
        d2xd2z: list = laplacian_didj(
            coord, 4, exp, center, lx, ly, lz, output, dalton_normalization, driver_time
        )
        d2yd2y: list = laplacian_didj(
            coord, 1, exp, center, lx, ly, lz, output, dalton_normalization, driver_time
        )
        d2yd2z: list = laplacian_didj(
            coord, 6, exp, center, lx, ly, lz, output, dalton_normalization, driver_time
        )
        d2zd2z: list = laplacian_didj(
            coord, 2, exp, center, lx, ly, lz, output, dalton_normalization, driver_time
        )
        cx: float = 3.0
        cy: float = 3.0
        cz: float = 3.0
        if component == 0:
            cx = 6.0
        elif component == 4:
            cy = 6.0
        else:
            cz = 6.0
        szmv: list = [
            cx * (a + b + c) + cy * (d + b + e) + cz * (c + e + f)
            for a, b, c, d, e, f in zip(d2xd2x, d2xd2y, d2xd2z, d2yd2y, d2yd2z, d2zd2z)
        ]
    elif component == 1 or component == 3:  # Lapdxdy or Lapdydx
        d2xdxdy = laplacian_didj(
            coord,
            15,
            exp,
            center,
            lx,
            ly,
            lz,
            output,
            dalton_normalization,
            driver_time,
        )
        d2ydxdy = laplacian_didj(
            coord,
            17,
            exp,
            center,
            lx,
            ly,
            lz,
            output,
            dalton_normalization,
            driver_time,
        )
        d2zdxdy = laplacian_didj(
            coord,
            14,
            exp,
            center,
            lx,
            ly,
            lz,
            output,
            dalton_normalization,
            driver_time,
        )
        szmv = [-(a + b + c) for a, b, c in zip(d2xdxdy, d2ydxdy, d2zdxdy)]
    elif component == 2 or component == 6:  # Lapdxdz or Lapdzdx
        d2xdxdz = laplacian_didj(
            coord,
            16,
            exp,
            center,
            lx,
            ly,
            lz,
            output,
            dalton_normalization,
            driver_time,
        )
        d2ydxdz = laplacian_didj(
            coord,
            11,
            exp,
            center,
            lx,
            ly,
            lz,
            output,
            dalton_normalization,
            driver_time,
        )
        d2zdxdz = laplacian_didj(
            coord,
            19,
            exp,
            center,
            lx,
            ly,
            lz,
            output,
            dalton_normalization,
            driver_time,
        )
        szmv = [-(a + b + c) for a, b, c in zip(d2xdxdz, d2ydxdz, d2zdxdz)]
    elif component == 5 or component == 7:  # dydz or dzdy
        d2xdydz = laplacian_didj(
            coord, 9, exp, center, lx, ly, lz, output, dalton_normalization, driver_time
        )
        d2ydydz = laplacian_didj(
            coord,
            18,
            exp,
            center,
            lx,
            ly,
            lz,
            output,
            dalton_normalization,
            driver_time,
        )
        d2zdydz = laplacian_didj(
            coord,
            20,
            exp,
            center,
            lx,
            ly,
            lz,
            output,
            dalton_normalization,
            driver_time,
        )
        szmv = [-(a + b + c) for a, b, c in zip(d2xdydz, d2ydydz, d2zdydz)]
    else:
        raise ValueError(
            f"***Error\n\n This term of Nabla . Nabla^T + Nabla^T . Nabla not exist: {component}"
        )

    if output > 10:
        driver_time.add_name_delta_time(
            name=f"Mass-Valocity correction to Spin--Zeeman {component}",
            delta_time=(time() - start),
        )

    return szmv


if __name__ == "__main__":
    # STO-2G
    print("\n LiH \n")
    s = szmv(
        coord=[[0.0, 0.0, -0.545857052], [0.0, 0.0, 2.309057052]],
        component=0,
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
        dalton_normalization=False,
        driver_time=None,
    )
    print("szmv : ", s, "\n", len(s), "\n\n")
