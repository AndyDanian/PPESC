import numpy as np

from molecule import *


class cluster():
    def __init__(
        self,
        coord: list[list[str]],
        basis: list[list[dict[str, list[float]]]],
    ):
        """
        Cluster object

        Args:
        -----
            coord (list[list]): atomic information of the different molecules
                [
                    ['element,charge,x,y,z', ...],
                    [...]
                ]
            basis (list[dict]): atomic basis set
                [
                    [
                        {'s':[exponents], 'p':[exponents], 'd':[exponents], 'f':[exponents]},
                        {'s':[exponents]}, {...}, ...
                    ],
                    [
                        ...
                    ]
                ]

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

        if len(np.array(coord).shape) > 2:
            raise ValueError(
                "*** Error \n\n\
                Systems with by above cluster doesn't exits.\n\
                Only exist:\n\
                        atoms (0D)<-molecule (1D)<-cluster (2D)\n\
                "
            )

        self._coord: list[list[str]] = coord
        self._basis: list[list[dict[str, list[float]]]] = basis

    ##################################################################
    # ATRIBUTES
    ##################################################################
    @property
    def molecules_number(self) -> int:
        "Molecule Number"
        if np.array(self._coord).ndim > 1:
            return len(self._coord)
        else:
            return 1

    @property
    def atoms_number(self) -> list:
        return [molecule(mol, self._basis[index]).atoms_number for index, mol in enumerate(self._coord)]

    @property
    def molecules_xyz(self) -> list:
        "Molecule Coordinates"
        return [molecule(coord, self._basis[index]).molecule_xyz for index, coord in enumerate(self._coord)]

    @property
    def atomic_symbols(self) -> list:
        "Atomic Symbols"
        return [molecule(coord, self._basis[index]).atomic_symbols for index, coord in enumerate(self._coord)]

    @property
    def atomic_numbers(self) -> list:
        return [molecule(coord, self._basis[index]).atomic_numbers for index, coord in enumerate(self._coord)]

    @property
    def charges(self) -> list:
        return [molecule(coord, self._basis[index]).charges for index, coord in enumerate(self._coord)]
    @property
    def primitives_number(self) -> int:
        return sum([molecule(coord, self._basis[index]).primitives_number for index, coord in enumerate(self._coord)])

    @property
    def angular_momentums(self) -> list:
        return [molecule(coord, self._basis[index]).angular_momentums for index, coord in enumerate(self._coord)]

    @property
    def exponents(self) -> list:
        return [molecule(coord, self._basis[index]).exponents for index, coord in enumerate(self._coord)]

    @property
    def mlx(self) -> list:
        return [molecule(coord, self._basis[index]).mlx for index, coord in enumerate(self._coord)]

    @property
    def mly(self) -> list:
        return [molecule(coord, self._basis[index]).mly for index, coord in enumerate(self._coord)]

    @property
    def mlz(self) -> list:
        return [molecule(coord, self._basis[index]).mlz for index, coord in enumerate(self._coord)]

    @property
    def amount_angular_momentums(self) -> list:
        "Amount of each Angular Momentum"
        cluster_angular_momentums = []
        for mol in self._basis:
            angular_momentums = []
            for basis in mol:
                l_momentums = {}
                for l, exp in basis.items():
                    l_momentums[l] = len(exp)
                angular_momentums.append(l_momentums)
            cluster_angular_momentums.append(angular_momentums)
        return cluster_angular_momentums

    ##################################################################
    # METHODS
    ##################################################################
    def get_atoms(self, verbose: int = 0):
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
        for index, mol in enumerate(self._coord):
            cluster_array.append(
                molecule(mol, self._basis[index]).get_atoms(verbose = verbose)
            )

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
        [
            {"s": [3.0, 0.01], "p": [1.0, 0.5]},
            {"s": [3.0, 0.01], "p": [1.0, 0.5]}
        ],
        [
            {"s": [3.0, 0.01], "p": [1.0, 0.5]},
            {"s": [3.0, 0.01], "p": [1.0, 0.5]},
        ]
        ]
    )

    print("\n # mol ", two_h2.molecules_number)
    print("# atoms ", two_h2.atoms_number)
    print("\ncluster   ", two_h2.molecules_xyz)
    print("\ncluster   ", two_h2.atomic_symbols)
    print("\ncluster   ", two_h2.atomic_numbers)
    print("\ncluster   ", two_h2.charges)
    print("\ncluster   ", two_h2.primitives_number)
    print("\ncluster   ", two_h2.angular_momentums)
    print("\ncluster   ", two_h2.exponents)
    print("\ncluster   ", two_h2.mlx)
    print("\ncluster   ", two_h2.mly)
    print("\ncluster   ", two_h2.mlz)
    print("\ncluster   ", two_h2.amount_angular_momentums)
    print("\ncluster   ", two_h2.get_atoms(verbose = 101))
