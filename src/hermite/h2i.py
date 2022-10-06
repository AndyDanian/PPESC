from libint import *


def h2i(
    name: str,
    charge: list[float] = [],
    coord: list[list[float]] = [],
    exp: list[float] = [],
    center: list[int] = [],
    lx: list[int] = [],
    ly: list[int] = [],
    lz: list[int] = [],
    verbose: int = 0,
    dalton_normalization: bool = None,
    driver_time: object = None,
):
    """
    Initialize the calculate of the integral choose

    Args:
        name (str): Integral name to calculate
        drive_time (drv_object): Object to manage the time
        charge (list): Atomic charge
        coord (list): Atomic Coordinate
        exp (list): Exponential of the primitives
        center (list): Where is centering of primitives
        lx, ly, lz (list): Component spatial of angular momentum component
        magnetic_xyz (int): Component spatial of the magnetic field
        spatial_sym (int): Spatial symmetry
        r_gauge (list): Gauge origen
        r_dipole (list): Dipole origen
        verbose (int): verbose level for integral calculation
    """
    if name.lower() == "e2pot":
        integral_twobody: list = e2pot(
            coord, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time
        )

    return integral_twobody
