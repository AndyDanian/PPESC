from lib import *
from h1i import *
from wave_function import *

nucleu = {"overlap": 0, "pot": 1}
esp_sym = {"overlap": 0, "pot": 0}
integral_symmetry = {"overlap": "sym", "pot": "sym"}


class eint:
    def __init__(self, wf: dict = None):
        """
        Manages the electronic integrals calculations

        Args:
            wf (dict): dictionary with information about wave function
        """

        if not wf:
            raise ValueError(
                "*** Error \n\n There isn't information in the  wave function object."
            )

        self._array_wf = wf
        # linealize arrays to calculate the integrals
        (
            self._charge,
            self._coord,
            self._exp,
            self._center,
            self._lx,
            self._ly,
            self._lz,
        ) = linealize_array_wf(wf)

    ##################################################################
    # METHODS
    ##################################################################

    def integration(
        self, atoms: list = None, names: list = None, output: int = None
    ):

        if not names:
            raise ValueError("***Error \n\n what integral do you want?")

        integrals = {}
        symmetry = {}

        for name in names:
            if nucleu[name.lower()] == 0 and esp_sym[name.lower()] == 0:
                symmetry[str(name)] = integral_symmetry[name.lower()]
                integrals[str(name)] = h1i(
                    self._charge,
                    self._coord,
                    self._exp,
                    self._center,
                    self._lx,
                    self._ly,
                    self._lz,
                    name,
                    output,
                    atoms,
                )
            if nucleu[name.lower()] == 1 and esp_sym[name.lower()] == 0:
                for atom in atoms:
                    symmetry[
                        name.lower()[0:3] + " " + str(atom + 1)
                    ] = integral_symmetry[name.lower()]
                    integrals[name.lower()[0:3] + " " + str(atom + 1)] = h1i(
                        self._charge,
                        self._coord,
                        self._exp,
                        self._center,
                        self._lx,
                        self._ly,
                        self._lz,
                        name,
                        output,
                        atom,
                    )

                    # Print integral
                    if output > 10:
                        print_triangle_matrix(
                            vector_to_matrix(
                                len(self._exp),
                                integrals[
                                    name.lower()[0:3] + " " + str(atom + 1)
                                ],
                                symmetry[
                                    name.lower()[0:3] + " " + str(atom + 1)
                                ],
                            ),
                            name.lower()[0:3] + " " + str(atom + 1),
                        )
        return integrals


if __name__ == "__main__":
    from wave_function import *

    wfn = wave_function("input/H2.molden")

    s = eint(wfn.build_wfn_array())

    s.integration([0, 1], ["pot"], 12)
