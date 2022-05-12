from atom import atom


class molecule(atom):
    def __init__(
        self,
        coord: list = None,
        basis: list = None,
    ):
        """
        Molecule object

        Args:

            coord (list[]): atomic information of a molecule
            coord = ['element,charge,x,y,z', 'element', '...']

            basis (list[dict]): basis set of atoms
            [{'s':[exponents], 'p':[exponents], 'd':[exponents], 'f':[exponents]},{'s':[exponents]}...]

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

        self._coord = coord
        self._basis = basis

        self._molecule_array = self.build_molecule_array(
            self._coord, self._basis
        )

    ##################################################################
    # METHODS
    ##################################################################

    def build_molecule_array(self, coord, basis):
        """
        Build one list of dictionaries with the atomic information

        [{"element": "symbol", "charge": 1.0, "x": 0.0, "y": 0.0, "z": 0.0, "mlx":[], "mly":[],
        "mlz":[], "exp":[]}, {"element": "symbol", ...}, ...]
        """

        molecule_array = []
        for count, atom in enumerate(coord):
            molecule_array.append(self.build_atom_array(atom, basis[count]))

        return molecule_array


if __name__ == "__main__":
    """
    Example to use molecule object
    """
    h2 = molecule(
        ["H 1.0 0.0 0.0 0.0", "H 1.0 0.0 0.0 0.75"],
        [
            {"s": [3.0, 0.01], "p": [1.0, 0.5]},
            {"s": [3.0, 0.01], "p": [1.0, 0.5]},
        ],
    )

    print("coordinate ", h2._coord)
    print("basis set ", h2._basis)
    print("Array ", h2._molecule_array)
