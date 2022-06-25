from atom import *

class molecule():
    def __init__(
        self,
        coord: list = None,
        basis: list = None,
    ):
        """
        Molecule object

        Args:
        ----
            coord (list[]): atomic information of a molecule
                            coord = ['element,charge,x,y,z', 'element', '...']
            basis (list[dict]): basis set of atoms
                                [{'s':[exponents], 'p':[exponents], ...},
                                    {'s':[exponents], ...}, ...]
        """

        if not coord:
            raise ValueError(
                "*** ERROR \n\n No input data \n\n\
            You can provide the coordinates, basis set, and molecular\
            orbitals information manually by mean of the list of \
            dictionaries coords, and basis.\n\
            --Minimal information is the coordinates."
            )
        # NOTE: The primitives are always taking as cartessians

        if np.array(coord).ndim > 1:
            raise ValueError(
                "*** Error \n\n\
                Systems with more one molecule is a cluster.\n\
                Use the cluster object:\n\
                        molecule/cluster(coord, basis)\n\
                "
            )

        self._coord = coord
        self._basis = basis
    ##################################################################
    # ATRIBUTES
    ##################################################################
    @property
    def atoms_number(self) -> int:
        "Atoms Number"
        return len(self._coord)

    @property
    def molecule_xyz(self) -> list:
        "Molecule Coordinates"
        return [atom(coord, self._basis[index]).atom_xyz for index, coord in enumerate(self._coord)]

    @property
    def atomic_symbols(self) -> list:
        "Atomic Symbols"
        return [atom(coord, self._basis[index]).atomic_symbol for index, coord in enumerate(self._coord)]

    @property
    def atomic_numbers(self) -> list:
        return [atom(coord, self._basis[index]).Z for index, coord in enumerate(self._coord)]

    @property
    def charges(self) -> list:
        return [atom(coord, self._basis[index]).q for index, coord in enumerate(self._coord)]
    @property
    def primitives_number(self) -> int:
        return sum([atom(coord, self._basis[index]).primitive_number for index, coord in enumerate(self._coord)])

    @property
    def angular_momentums(self) -> list:
        return [atom(coord, self._basis[index]).angular_momentum for index, coord in enumerate(self._coord)]

    @property
    def exponents(self) -> list:
        return [atom(coord, self._basis[index]).exponents for index, coord in enumerate(self._coord)]

    @property
    def mlx(self) -> list:
        return [atom(coord, self._basis[index]).mlx for index, coord in enumerate(self._coord)]

    @property
    def mly(self) -> list:
        return [atom(coord, self._basis[index]).mly for index, coord in enumerate(self._coord)]

    @property
    def mlz(self) -> list:
        return [atom(coord, self._basis[index]).mlz for index, coord in enumerate(self._coord)]

    @property
    def amount_angular_momentums(self) -> list:
        "Amount of each Angular Momentum"
        angular_momentums = []
        for basis in self._basis:
            l_momentums = {}
            for l, exp in basis.items():
                l_momentums[l] = len(exp)
            angular_momentums.append(l_momentums)
        return angular_momentums

    ##################################################################
    # METHODS
    ##################################################################
    #def build_molecule_array(self, verbose: int = None):
    def get_atoms(self, verbose: int = None):
        """
        Build one list of dictionaries with the atomic information

        [
            {
                "element": "symbol", "charge": 1.0, "x": 0.0, "y": 0.0, "z": 0.0,
                "mlx":[...], "mly":[...], "mlz":[...], "exp":[]
            },
            {
                "element": "symbol", ...
            },
            ...
        ]
        """

        molecule_array: list = []
        for count, atomic_information in enumerate(self._coord):
            molecule_array.append(atom(atomic_information, self._basis[count]).build_atom_array(verbose=verbose))

        return molecule_array


if __name__ == "__main__":
    """
    Example to use molecule object
    """
    h2 = molecule(
        ["He 1.0 0.0 0.0 0.0", "H 1.0 0.0 0.0 0.75"],
        [
            {"s": [3.0, 0.01], "p": [1.0, 0.5]},
            {"s": [3.0, 0.01], "p": [1.0, 0.5]},
        ],
    )

    h2.get_atoms(verbose=101)
