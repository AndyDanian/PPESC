from molecule import molecule


class cluster(molecule):
    def __init__(
        self,
        coord: list = None,
        basis: list = None,
    ):
        """
        Cluster object

        Args:

            coord (list[list]): atomic information of the different molecules
            coord = [['element,charge,x,y,z', ...], [...],...]

            basis (list[dict]): atomic basis set
            [{'s':[exponents], 'p':[exponents], 'd':[exponents], 'f':[exponents]},
            {'s':[exponents]}, {...}, ...]

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

        self._cluster_array = self.build_cluster_array(
            self._coord, self._basis
        )

    def build_cluster_array(self, coord, basis):
        """
        Build one list of lists of dictionaries with the molecule information

        [
            [
                {"element": "symbol", "charge": 1.0, "x": 0.0, "y": 0.0, "z": 0.0,
                "mlx":[], "mly":[], "mlz":[], "exp":[]},
                {"element": "symbol", ... }...
            ],
            [
                {"element": "symbol", ...}
            ], ...
        ]
        """

        cluster_array = []
        init_atom = 0
        end_atom = 0
        for molecule in coord:
            end_atom += len(molecule)
            cluster_array.append(
                self.build_molecule_array(molecule, basis[init_atom:end_atom])
            )
            init_atom = end_atom

        return cluster_array


if __name__ == "__main__":
    """
    Example to use cluster object
    """
    two_h2 = cluster(
        [
            ["H 1.0 0.0 0.0 0.0", "H 1.0 0.0 0.0 0.75"],
            ["H 1.0 1.0 0.0 0.0", "H 1.0 1.0 0.0 0.75"],
        ],
        [
            {"s": [3.0, 0.01], "p": [1.0, 0.5]},
            {"s": [3.0, 0.01], "p": [1.0, 0.5]},
            {"s": [3.0, 0.01], "p": [1.0, 0.5]},
            {"s": [3.0, 0.01], "p": [1.0, 0.5]},
        ],
    )

    print("\ncoordinate ", two_h2._coord)
    print("basis set ", two_h2._basis)
    print("\ncluster   ", two_h2._cluster_array)
