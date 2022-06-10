from libsrc import * # from io import molden as mn
from cluster import *


def linealize_array_wf(wfn_array: dict or list = None):
    """
    Linealize the array into wave function like exponentials
    to calculate the integrals

    Return:
        charge (list): list 1d of charges
        coord (list): list 2d with coordinates of the atoms
        exp (list): list 1d with the exponentials
        center (list): list 1d with the center of the gaussian
        lx (list): list 1d with the x component of ml of the gaussian
        ly (list): list 1d with the y component of ml of the gaussian
        lz (list): list 1d with the z component of ml of the gaussian
    """

    if isinstance(wfn_array, dict):
        cartessian_primitives = wfn_array["cto"]
        wfn_array = wfn_array["cluster"]

    charge = [q["charge"] for atom in wfn_array for q in atom]
    coord = [[r["x"], r["y"], r["z"]] for atom in wfn_array for r in atom]
    exp = [i for atom in wfn_array for e in atom for i in e["exp"]]
    center = [
        count
        for count, atom in enumerate(wfn_array)
        for e in atom
        for i in e["exp"]
    ]
    mlx = [mlx for atom in wfn_array for ml in atom for mlx in ml["mlx"]]
    mly = [mly for atom in wfn_array for ml in atom for mly in ml["mly"]]
    mlz = [mlz for atom in wfn_array for ml in atom for mlz in ml["mlz"]]
    angular_moments = [symbol for atom in wfn_array for l in atom for symbol, amount in l["l"].items() for i in range(amount)]

    return charge, coord, exp, center, mlx, mly, mlz, angular_moments, cartessian_primitives


class wave_function(cluster):
    def __init__(
        self,
        filename: str = None,
        coord: list = None,
        basis: list = None,
        mos: list = None,
        cartessian_primitive: bool = False,
    ):
        """
        Wave function object need a filename with wave function

        Args:
            filename (string): file with wave function (molden/wfnx)

            coord (list[list]): atomic information of the different molecules
            coord = [['element,charge,x,y,z', ...], [...],...]

            basis (list[dict]): atomic basis set
            [{'s':[exponents], 'p':[exponents], 'd':[exponents], 'f':[exponents]},{'s':[exponents]}...]

            mos (list[dict]): molecular orbitals information (coefficients, ...)
            [{'energy':..., 'spin':..., "occupation":..., "coefficients":[...]},...]

            cartessian_primitive (bool): if True, the molecular orbitals are in cartessian
        """

        if not filename:
            if not coord or not basis or not mos or not cartessian_primitive:
                raise ValueError(
                    "*** ERROR \n\n No input data \n\n\
                You can used a molden or wfnx file. Also, you can provide the\
                coordinates, basis set, and molecular orbitals information\
                manually by mean of the list of dictionaries coords, basis,\
                mos, and cartessian_primitive."
                )

        if filename:
            # get information from molden/wfnx file
            (
                self.coord,
                self.basis,
                self.mos,
                self.cartessian_primitive,
            ) = read_molden(
                filename,
            )
        else:
            self.coord = coord
            self.basis = basis
            self.mos = mos
            self.cartessian_primitive = cartessian_primitive

    ##################################################################
    # METHODS
    ##################################################################

    def build_wfn_array(self):
        """
        Build one list of lists of dictionaries with the molecule information and MOs.
        This method return a wave function encapsuled, i.e., all information in one array

        {
            "cluster":[
                    [
                        {"element": "symbol", "charge": 1.0, "x": 0.0, "y": 0.0, "z": 0.0,
                        "l": {"s": 5, "p": 4, ...}, "mlx":[], "mly":[], "mlz":[], "exp":[]},
                        {"element": "symbol", ...}], ...
                    ],
            "mos":[
                    {'energy':..., 'spin':...,"occupation":..., 'coefficients':[...]},
                    {'energy':..., ...}, ...
                ], ...
        }

        Cluster is an object with minimun one molecule or one atom
        """

        wfn = {}
        molecule = []

        molecule = [
            x for x in self.build_cluster_array(self.coord, self.basis)
        ]

        wfn["cluster"] = molecule
        wfn["mos"] = self.mos
        wfn["cto"] = self.cartessian_primitive

        return wfn


if __name__ == "__main__":
    """
    Example to use wave function object
    """
    wfn = wave_function("io/H2O.molden")

    array = wfn.build_wfn_array()
    print("\nWave function ", array["cluster"][0])
    print("\nmlx amount ", len(array["cluster"][0][0]["mlx"]))
    print("\nmly amount ", len(array["cluster"][0][0]["mly"]))
    print("\nmlz amount ", len(array["cluster"][0][0]["mlz"]))
    print("\nexp amount ", len(array["cluster"][0][0]["exp"]))

    print("\n Wafe Function like dictionary ")
    charge, coord, exp, center, lx, ly, lz = linealize_array_wf(array)
    print(" charge ", charge)

    print("\n Wafe Function like list ")
    charge, coord, exp, center, lx, ly, lz = linealize_array_wf(
        array["cluster"]
    )
    print(" charge ", charge)
