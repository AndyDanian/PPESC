from libsrc import *

class atom:
    def __init__(self, coord: str or list = None, basis: dict = None):
        """
        Atom object

        Args:
            coord (string, list): atomic information
            coord = 'H charge x y z'

            basis (dict): keywords are the angular quatum number and values are
            atomic exponents basis set
            {'s':[exponents], 'p':[exponents], 'd':[exponents], 'f':[exponents]}
        """

        if not coord or not basis:
            raise ValueError(
                "*** ERROR \n\n No atomic data \n\n\
                You must provide the coordinates and exponents of basis set.\n\
                Example: atom('H 1.0 0.0 0.0 0.0', '{'s:'3.0,2.0''}')"
            )

        self._coord = coord
        self._basis = basis
        self._atom_array = self.build_atom_array(self._coord, self._basis)

    ##################################################################
    # METHODS
    ####################################basis##############################
    def build_atom_array(self, coord, basis):
        """
        Build a dictionary for atom

        {"element":"H", "charge":1.0, "x":0.0, "y":0.0, "z":0.0, l:{"s":5, "p":4, ..},
        "mlx":[], "mly":[], "mlz":[], "exp":[] }
        """

        if isinstance(coord, list):
            if len(coord) > 1:
                raise ValueError(
                    "*** Error \n\
                        In the atom object, only accept list with one element"
                )
            coord = coord[0]

        atom_array = {}
        atom_array["element"] = coord.split()[0]

        atom_array["charge"] = float(coord.split()[1])

        atom_array["x"] = float(coord.split()[2])
        atom_array["y"] = float(coord.split()[3])
        atom_array["z"] = float(coord.split()[4])

        angular_q = {}
        mlx = []
        mly = []
        mlz = []
        exponents = []
        for l, exp in basis.items():
            angular_q[l] = len(exp)
            mlx += cartessian_mlx[l] * len(exp)
            mly += cartessian_mly[l] * len(exp)
            mlz += cartessian_mlz[l] * len(exp)
            exponents += exp * angular_number[l]
        atom_array["l"] = angular_q
        atom_array["mlx"] = mlx
        atom_array["mly"] = mly
        atom_array["mlz"] = mlz
        atom_array["exp"] = exponents

        return atom_array


if __name__ == "__main__":
    """
    Example to use atom object

    S: Atomic symbol
    Z: Atomic number or charge
    l: angular quantum number
    exp: Exponent value for gaussian
    """
    #               [ S  Z   X   Y   Z  ]    l    exp         l     exp
    hydrogen = atom(["H 1.0 0.0 0.0 0.0"], {"s": [3.0, 0.01], "p": [1.0, 0.5]})

    print(" Coordinate is a list ")
    print("\ncoordinate ", hydrogen._coord)
    print("basis set ", hydrogen._basis)
    print("Array ", hydrogen._atom_array)

    hydrogen = atom("H 1.0 0.0 0.0 0.0", {"s": [3.0, 0.01], "p": [1.0, 0.5]})
    print(" Coordinate is a string ")
    print("\ncoordinate ", hydrogen._coord)
    print("basis set ", hydrogen._basis)
    print("Array ", hydrogen._atom_array)
