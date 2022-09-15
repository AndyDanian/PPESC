from libint import *
from h2int.e2pot import *

def h2i(
    charge: list = None,
    coord: list = None,
    exp: list = None,
    center: list = None,
    lx: list = None,
    ly: list = None,
    lz: list = None,
    name: str = None,
    verbose: int = 0,
    dalton_normalization: bool = None,
    driver_time: object = None
):

    if name.lower() == "e2pot":
        integral_twobody: list = e2pot(coord, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time)

    return integral_twobody