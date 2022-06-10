from libint import *


def h2i(
    charge: list = None,
    coord: list = None,
    exp: list = None,
    center: list = None,
    lx: list = None,
    ly: list = None,
    lz: list = None,
    name: str = None,
    output: int = 0,
    dalton_normalization: bool = None,
):

    if name.lower() == "e2pot":
        integral_twobody: list = e2pot(coord, exp, center, lx, ly, lz, output, dalton_normalization)

    return integral_twobody