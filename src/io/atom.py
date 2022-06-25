from libsrc import *

class atom:
    def __init__(self, coord: str or list = None, basis: dict = None):
        """
        Atom object

        Args:
        ----
            coord (string, list): atomic information
                                    coord = 'H charge x y z'
            basis (dict): keywords are the angular quantum number and values are
                            atomic exponents basis set
                            {'s':[exponents], 'p':[exponents], 'd':[exponents], ...}
        """

        if not coord or not basis:
            raise ValueError(
                "*** ERROR \n\n No atomic data \n\n\
                You must provide the coordinates and exponents of basis set.\n\
                Example: atom('H 1.0 0.0 0.0 0.0', '{'s:'3.0,2.0''}')"
            )

        if np.array(coord).ndim > 1 or (isinstance(coord, list) and len(coord) > 1):
            raise ValueError(
                "*** Error \n\n\
                Systems with more one atom/molecule is a cluster.\n\
                Use the molecule/cluster object:\n\
                        molecule/cluster(coord, basis)\n\
                "
            )

        if (isinstance(coord, str) and len(coord.split(" ")) != 5) or (isinstance(coord, list) and len(coord[0].split(" ")) != 5):
            raise ValueError(
            "*** Error \n\n\
            The coordinate vairable is conformed by 5 strings: label or atomic number,\n\
            atomic charge, and cartessain coordinates\n. \
                \"H 1.0 0.00 0.00 0.00\"\n\
            "
            )

        self._coord = coord
        self._basis = basis

    ##################################################################
    # ATTRIBUTES
    ##################################################################
    @property
    def atom_xyz(self) -> list:
        "Return Coordinates"
        if isinstance(self._coord, list):
            coordinates = self._coord[0].split(" ")
        elif isinstance(self._coord, str):
            coordinates = self._coord.split(" ")
        return [float(r) for i, r in enumerate(coordinates) if i > 1]

    @property
    def atomic_symbol(self) -> str:
        "Atom Symbol"
        if isinstance(self._coord, list):
            if self._coord[0].split(" ")[0].isalpha():
                return self._coord[0].split(" ")[0]
            elif self._coord[0].split(" ")[0].isnumeric():
                return atomic_symbol[self._coord[0].split(" ")[0]]
        elif isinstance(self._coord, str):
            if self._coord.split(" ")[0].isalpha():
                return self._coord.split(" ")[0]
            elif self._coord.split(" ")[0].isnumeric():
                return atomic_symbol[self._coord[0].split(" ")[0]]

    @property
    def atomic_number(self) -> int:
        "Atomic Number"
        if isinstance(self._coord, list):
            if self._coord[0].split(" ")[0].isalpha():
                return [Z for Z, symbol in atomic_symbol.items() if symbol == self._coord[0].split(" ")[0]][0]
            elif self._coord[0].split(" ")[0].isnumeric():
                return int(self._coord[0].split(" ")[0])
        elif isinstance(self._coord, str):
            if self._coord.split(" ")[0].isalpha():
                return [int(Z) for Z, symbol in atomic_symbol.items() if symbol == self._coord.split(" ")[0]][0]
            elif self._coord.split(" ")[0].isnumeric():
                return int(self._coord.split(" ")[0])

    @property
    def Z(self) -> int:
        return self.atomic_number

    @property
    def charge(self) -> float:
        "Atomic Charge"
        if isinstance(self._coord, list):
            q = self._coord[0].split(" ")[1]
        elif isinstance(self._coord, str):
            q = self._coord.split(" ")[1]
        return float(q)

    @property
    def q(self) -> float:
        return self.charge

    @property
    def angular_momentum(self) -> list:
        "Atomic Primitive Type"
        return [prim_type for prim_type in self._basis.keys()]

    @property
    def amount_angular_momentum(self) -> dict:
        "Amount of each Angular Momentum"
        angular_momentums = {}
        for l, exp in self._basis.items():
            angular_momentums[l] = len(exp)
        return angular_momentums

    @property
    def primitive_number(self) -> int:
        "Primitive Number"
        return sum([len(exp)*angular_number[l] for l, exp in self._basis.items()])

    @property
    def exponents(self) -> list:
        "Exponents"
        return [exp for l, exps in self._basis.items() for exp in exps for i in range(angular_number[l])]

    @property
    def mlx(self) -> list:
        "Exponent in the X direction"
        return [mlx for l, exps in self._basis.items() for mlx in cartessian_mlx[l] * len(exps)]

    @property
    def mly(self) -> list:
        "Exponent in the X direction"
        return [mly for l, exps in self._basis.items() for mly in cartessian_mly[l] * len(exps)]

    @property
    def mlz(self) -> list:
        "Exponent in the X direction"
        return [mlz for l, exps in self._basis.items() for mlz in cartessian_mlz[l] * len(exps)]

    ##################################################################
    # METHODS
    ##################################################################
    def build_atom_array(self, verbose: int = None):
        """
        Build a dictionary for atom

        {"element":"H", "charge":1.0, "x":0.0, "y":0.0, "z":0.0, l:{"s":5, "p":4, ..},
        "mlx":[], "mly":[], "mlz":[], "exp":[] }
        """

        atom_array = {}
        atom_array["element"] = self.atomic_symbol

        atom_array["charge"] = float(self.q)

        atom_array["x"] = float(self.atom_xyz[0])
        atom_array["y"] = float(self.atom_xyz[1])
        atom_array["z"] = float(self.atom_xyz[2])

        atom_array["l"] = self.amount_angular_momentum
        atom_array["mlx"] = self.mlx
        atom_array["mly"] = self.mly
        atom_array["mlz"] = self.mlz
        atom_array["exp"] = self.exponents

        if isinstance(verbose, int) and verbose >= 100:
            print("\nAtomic Infomarion : ",atom_array["element"])
            print("*) Atomic number ",self.Z)
            print("*) Charge ",self.q)
            print("*) Coordinate ",self.atom_xyz)
            print("*) Angular Momentums ",self.amount_angular_momentum)
            print("*) Exponents ",self.exponents)
            print("*) mlx ",self.mlx)
            print("*) mly ",self.mly)
            print("*) mlz ",self.mlz)

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
    hydrogen = atom(["H 1.0 0.0 0.0 0.0"], {"s": [3.0, 0.01], "p": [1.0, 0.5, 2.0]})
    hydrogen = atom("1 1.0 0.0 0.0 0.0", {"s": [3.0, 0.01], "p": [1.0, 0.5, 2.0]})

    hydrogen.build_atom_array(verbose=101)