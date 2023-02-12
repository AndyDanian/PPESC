from typing import Union

from libsrc import *  # from io import molden as mn


class wave_function:
    def __init__(
        self,
        filename: str,
        scratch_path: Union[Path, str],
        job_folder: str = "",
        coord: list[list[str]] = [],
        basis: list[list[dict[str, list[float]]]] = None,
        mos: list[dict] = [],
        cartessian_primitive: bool = False,
        restart: int = 0,
    ):
        """
        Wave function object need a filename with wave function

        Args:
            filename (string): file with wave function (molden/wfnx)
            coord (list[list]): atomic information of the different molecules
                coord = [
                            ['element,charge,x,y,z', ...],
                            [...],...
                        ]
            basis (list[dict]): atomic basis set
                basis = [
                            [
                                {'s':[exponents], 'p':[exponents], 'd':[exponents], 'f':[exponents]},
                            ...]
                            [
                                {'s':[exponents]},...
                            ]
                        ]

            mos (list[dict]): molecular orbitals information (coefficients, ...)
            [{'energy':..., 'spin':..., "occupation":..., "coefficients":[...]},...]

            cartessian_primitive (bool): if True, the molecular orbitals are in cartessian
            restart (int): Errase (0) or not (1) information into de schratch directory (default: 0)
        """

        start = time()
        self._driver_time: drv_time = drv_time()
        self._driver_scratch: scratch = scratch(
            scratch=scratch_path, job_folder=job_folder, restart=restart
        )

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
            file_output: str = (
                str(self._driver_scratch.job_path)
                + "/"
                + Path(filename).stem.upper()
                + ".out"
            )
            self._driver_scratch.output_path = Path(file_output)
            self._driver_scratch.write_header_output()
            # get information from molden/wfnx file
            (
                self._coord,
                self._basis,
                self._mos,
                self._cartessian_primitive,
            ) = read_molden(file_molden=filename, drv_scratch=self._driver_scratch)
        else:
            self._driver_scratch._output_path = Path("OUTPUT.out")
            self._driver_scratch.write_header_output()
            self._coord = coord
            self._basis = basis
            self._mos = mos
            self._cartessian_primitive = cartessian_primitive

        self.wave_function_information()
        self._driver_time.add_name_delta_time(
            name="Reading wave function", delta_time=(time() - start)
        )

        self._driver_scratch.write_output(type=2, drv_time=self._driver_time)
        self._driver_time.reset
        self._driver_scratch.write_output("\n")

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
    def coordinates(self) -> list:
        "Atoms coordinates"
        return [
            [float(at.split()[2]), float(at.split()[3]), float(at.split()[4])]
            for mol in self._coord
            for at in mol
        ]

    @property
    def atom_number(self) -> int:
        "Atom Number"
        return len(self.coordinates)

    @property
    def atomic_numbers(self) -> Union[None, list]:
        "Atomic Numbers"
        Z: list = []
        for mol in self._coord:
            for at in mol:
                if at.split(" ")[0].isalpha():
                    Z.append(
                        [
                            int(iZ)
                            for iZ, symbol in atomic_symbol.items()
                            if symbol == at.split(" ")[0]
                        ][0]
                    )
                elif self._coord.split()[0].isnumeric():
                    return Z.append(int(self._coord.split()[0]))
        return Z

    @property
    def mo_occ(self) -> int:
        "Molecular Orbital Occupied Number"
        return len([mo["occupation"] for mo in self._mos if mo["occupation"] > 0])

    @property
    def charges(self) -> list[float]:
        "Atomic Charges"
        return [float(at.split()[1]) for mol in self._coord for at in mol]

    @property
    def atomic_symbols(self) -> list:
        "Atomic Symbols"
        return [
            at.split()[0]
            if isinstance(at.split()[0], str)
            else atomic_symbol[int(at.split()[0])]
            for mol in self._coord
            for at in mol
        ]

    @property
    def cto(self) -> bool:
        "Primitive Symmetry"
        return self._cartessian_primitive

    @property
    def exponents(self) -> list:
        "Exponents"
        return [
            exp
            for mol in self._basis
            for at in mol
            for l, exps in at.items()
            for exp in exps
            for i in range(angular_number[l])
        ]

    @property
    def primitives_number(self) -> int:
        "Primitive number"
        return len(self.exponents)

    @property
    def primitives_number_sph(self) -> int:
        "Primitive number"
        return len(
            [
                exp
                for mol in self._basis
                for at in mol
                for l, exps in at.items()
                for exp in exps
                for i in range(angular_number_sph[l])
            ]
        )

    @property
    def primitives_number_car(self) -> int:
        "Primitive number"
        return len(
            [
                exp
                for mol in self._basis
                for at in mol
                for l, exps in at.items()
                for exp in exps
                for i in range(angular_number[l])
            ]
        )

    @property
    def mo_virt(self) -> int:
        "Molecular Orbital Occupied Number"
        return len([mo["occupation"] for mo in self._mos if mo["occupation"] == 0.0])

    @property
    def primitives_centers(self) -> list:
        "Primitive Center"
        center = []
        count = 0
        for mol in self._basis:
            for at in mol:
                for l, exps in at.items():
                    center += [count] * angular_number[l] * len(exps)
            count += 1
        return center

    @property
    def mlx(self) -> list:
        "Exponent in the X direction"
        return [
            mlx
            for mol in self._basis
            for at in mol
            for l, exps in at.items()
            for mlx in cartessian_mlx[l] * len(exps)
        ]

    @property
    def mly(self) -> list:
        "Exponent in the X direction"
        return [
            mly
            for mol in self._basis
            for at in mol
            for l, exps in at.items()
            for mly in cartessian_mly[l] * len(exps)
        ]

    @property
    def mlz(self) -> list:
        "Exponent in the X direction"
        return [
            mlz
            for mol in self._basis
            for at in mol
            for l, exps in at.items()
            for mlz in cartessian_mlz[l] * len(exps)
        ]

    @property
    def angular_momentums(self) -> list:
        "Atomic Primitive Type"
        return [
            l
            for mol in self._basis
            for at in mol
            for l, exps in at.items()
            for exp in exps
        ]

    @property
    def amount_angular_momentums(self) -> dict:
        "Amount of each Angular Momentum"
        angular_momentums: dict = {}
        for mol in self._basis:
            for at in mol:
                for l, exp in at.items():
                    if l in angular_momentums.keys():
                        angular_momentums[l] += len(exp)
                    else:
                        angular_momentums[l] = len(exp)
        return angular_momentums

    @property
    def amount_angular_momentums_by_atom(self) -> list:
        "Amount of each Angular Momentum"
        l_atom: list = []
        for mol in self._basis:
            for at in mol:
                angular_momentums: dict = {}
                for l, exp in at.items():
                    if l in angular_momentums.keys():
                        angular_momentums[l] += len(exp)
                    else:
                        angular_momentums[l] = len(exp)
                l_atom.append(angular_momentums)
        return l_atom

    @property
    def mo_coefficients(self) -> list:
        "Molecular Orbital Coefficients"
        return [mo["coefficients"] for mo in self._mos]

    @property
    def mo_energies(self) -> list:
        "Molecular Orbitals Energies"
        return [mo["energy"] for mo in self._mos]

    @mo_energies.setter
    def mo_energies(self, moe: list) -> list:
        "Molecular Orbitals Energies"
        self._mos = [{name: (moe[count]
                      if name == "energy"
                      else value)
                      for name, value in mo.items()
                      }
                      for count, mo in enumerate(self._mos)
                    ]
    ##################################################################
    # METHODS
    ##################################################################
    def wave_function_information(self):
        """
        Print information about wave function and molecule
        """

        self._driver_scratch.write_output("System ")
        self._driver_scratch.write_output("-" * 50)
        if np.max(self.coordinates) > 9999:
            form: str = "{:.4e}"
        else:
            form: str = "{:.4f}"
        len_s = [len(s) for s in self.atomic_symbols]
        if np.max(len_s) > 5:
            forms: str = "{:" + str(np.max(len_s)) + "s}"
        else:
            forms: str = "{:5s}"
        for s, xyz in zip(self.atomic_symbols, self.coordinates):
            information: str = "        " + forms.format(s).center(5)
            for x in xyz:
                information += str(form.format(x)).center(9)
            self._driver_scratch.write_output(information)
        self._driver_scratch.write_output("-" * 50 + "\n")

        self._driver_scratch.write_output("Primitive Informaiton")
        self._driver_scratch.write_output("-" * 50)
        total = 0
        count = 0
        for s, b in zip(self.atomic_symbols, self.amount_angular_momentums_by_atom):
            information: str = s.center(7) + str(self.charges[count]).center(5) + "  "
            count += 1
            for n, c in b.items():
                information += str(c) + n
            self._driver_scratch.write_output(information)
            total += sum(b.values())
        self._driver_scratch.write_output("=" * 50)
        self._driver_scratch.write_output(
            "Total: " + str(sum(self.charges)).center(5) + "   " + str(total) + "\n"
        )

        sample: str = "Cartessian Primitive: "
        self._driver_scratch.write_output(sample + str(self.primitives_number_car))
        self._driver_scratch.write_output(
            "Spherical Primitive: ".ljust(len(sample)) + str(self.primitives_number_sph)
        )
        self._driver_scratch.write_output(
            "Occupied Orbitals: ".ljust(len(sample)) + str(self.mo_occ)
        )
        self._driver_scratch.write_output(
            "Virtuals Orbitals: ".ljust(len(sample)) + str(self.mo_virt)
        )
        self._driver_scratch.write_output(
            "Ocuppied ⇌ Virtuals: ".ljust(len(sample)) + str(self.mo_occ * self.mo_virt)
        )
        self._driver_scratch.write_output("-" * 50 + "\n")


if __name__ == "__main__":
    """
    Example to use wave function object
    """
    wfn = wave_function("../tests/molden_file/H2.molden", scratch_path="/home1/scratch")

    print(" Molecule Number ", wfn.molecules_number)
    print(" Atom Number ", wfn.atom_number)
    print(" MO Occupied/Virtuals ", wfn.mo_occ, wfn.mo_virt)
    print(" Atomic Symbols ", wfn.atomic_symbols)
    print(" Coordinates ", wfn.coordinates)
    print(" Z ", wfn.atomic_numbers)
    print(" Charges ", wfn.charges)
    print(" Cartessian Symmetry ", wfn.cto)
    print(" Primitive Number ", wfn.primitives_number)
    print(" Angular Momentums ", wfn.angular_momentums)
    print(" Amount Angular Momentums ", wfn.amount_angular_momentums)
    print(" Exponents ", wfn.exponents)
    print(" Primitive Center ", wfn.primitives_centers)
    print(" mlx ", wfn.mlx)
    print(" mly ", wfn.mly)
    print(" mlz ", wfn.mlz)

    print("mos ",[{name: (value if name != 'energy' else 9999999.0) for name, value in mo.items()} for mo in wfn._mos])

    wfn._driver_scratch.remove_job_folder()
