from libsrc import * # from io import molden as mn

class wave_function():
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
        """

        print("*"*80)
        print("PROGRAM NAME".center(80))
        print("version 0.0".center(80))
        print("2022".center(80))
        print("Authors:".ljust(80))
        print("         Andy Zapata".ljust(80))
        print("*"*80)
        print()

        start = time()
        self._driver_time: drv_time = drv_time()

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
                self._coord,
                self._basis,
                self._mos,
                self._cartessian_primitive,
            ) = read_molden(
                filename,
            )
        else:
            self._coord = coord
            self._basis = basis
            self._mos = mos
            self._cartessian_primitive = cartessian_primitive

        self.wave_function_information()
        print_time(name = "Reading wave function", delta_time = time() - start, header = False, tailer =False)
        print()

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
        return [[float(at.split()[2]), float(at.split()[3]), float(at.split()[4])] for mol in self._coord for at in mol]

    @property
    def atom_number(self) -> int:
        "Atom Number"
        return len(self.coordinates)

    @property
    def atomic_numbers(self) -> list:
        "Atomic Numbers"
        Z = []
        for mol in self._coord:
            for at in mol:
                if at.split(" ")[0].isalpha():
                    Z.append([int(iZ) for iZ, symbol in atomic_symbol.items() if symbol == at.split(" ")[0]][0])
                elif self._coord.split()[0].isnumeric():
                    return Z.append(int(self._coord.split()[0]))
        return Z

    @property
    def mo_occ(self) -> int:
        "Molecular Orbital Occupied Number"
        return len([mo["occupation"] for mo in self._mos if mo["occupation"] > 0])

    @property
    def charges(self) -> list:
        "Atomic Charges"
        return [int(at.split()[1]) for mol in self._coord for at in mol]

    @property
    def atomic_symbols(self) -> list:
        "Atomic Symbols"
        return [at.split()[0] if isinstance(at.split()[0], str) else atomic_symbol(int(at.split()[0])) for mol in self._coord for at in mol]

    @property
    def cto(self) -> bool:
        "Primitive Symmetry"
        return self._cartessian_primitive

    @property
    def exponents(self) -> list:
        "Exponents"
        return [exp for mol in self._basis for at in mol for l, exps in at.items() for exp in exps for i in range(angular_number[l])]

    @property
    def primitives_number(self) -> int:
        "Primitive number"
        return len(self.exponents)

    @property
    def primitives_number_sph(self) -> int:
        "Primitive number"
        return len([exp for mol in self._basis for at in mol for l, exps in at.items() for exp in exps for i in range(angular_number_sph[l])])


    @property
    def primitives_number_car(self) -> int:
        "Primitive number"
        return len([exp for mol in self._basis for at in mol for l, exps in at.items() for exp in exps for i in range(angular_number[l])])

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
                    center += [count]*angular_number[l]*len(exps)
            count += 1
        return center

    @property
    def mlx(self) -> list:
        "Exponent in the X direction"
        return [mlx for mol in self._basis for at in mol for l, exps in at.items() for mlx in cartessian_mlx[l] * len(exps)]

    @property
    def mly(self) -> list:
        "Exponent in the X direction"
        return [mly for mol in self._basis for at in mol for l, exps in at.items() for mly in cartessian_mly[l] * len(exps)]

    @property
    def mlz(self) -> list:
        "Exponent in the X direction"
        return [mlz for mol in self._basis for at in mol for l, exps in at.items() for mlz in cartessian_mlz[l] * len(exps)]

    @property
    def angular_momentums(self) -> list:
        "Atomic Primitive Type"
        return [l for mol in self._basis for at in mol for l, exps in at.items() for exp in exps]

    @property
    def amount_angular_momentums(self) -> dict:
        "Amount of each Angular Momentum"
        angular_momentums = {}
        for mol in self._basis:
            for at in mol:
                for l, exp in at.items():
                    if l in angular_momentums.keys():
                        angular_momentums[l] += len(exp)
                    else:
                        angular_momentums[l] = len(exp)
        return angular_momentums

    @property
    def amount_angular_momentums_by_atom(self) -> dict:
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
    ##################################################################
    # METHODS
    ##################################################################
    def wave_function_information(self):
        """
        Print information about wave function and molecule
        """

        print("System ")
        print("-"*50)
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
            print(f"        ", forms.format(s).center(5), *[str(form.format(x)).center(9) for x in xyz])
        print("-"*50)
        print()

        print("Primitive Informaiton")
        print("-"*50)
        total = 0
        for s, b in zip(self.atomic_symbols, self.amount_angular_momentums_by_atom):
            print(forms.format(s).center(5), *[str(c) + n for n, c in b.items()])
            total += sum(b.values())
        print("="*50)
        print("Total: ",total)
        print("-"*50)
        print()

if __name__ == "__main__":
    """
    Example to use wave function object
    """
    wfn = wave_function("../tests/molden_file/H2.molden")

    print(" Molecule Number ",wfn.molecules_number)
    print(" Atom Number ",wfn.atom_number)
    print(" MO Occupied/Virtuals ",wfn.mo_occ,wfn.mo_virt)
    print(" Atomic Symbols ",wfn.atomic_symbols)
    print(" Coordinates ",wfn.coordinates)
    print(" Z ",wfn.atomic_numbers)
    print(" Charges ",wfn.charges)
    print(" Cartessian Symmetry ",wfn.cto)
    print(" Primitive Number ",wfn.primitives_number)
    print(" Angular Momentums ",wfn.angular_momentums)
    print(" Amount Angular Momentums ",wfn.amount_angular_momentums)
    print(" Exponents ",wfn.exponents)
    print(" Primitive Center ",wfn.primitives_centers)
    print(" mlx ",wfn.mlx)
    print(" mly ",wfn.mly)
    print(" mlz ",wfn.mlz)
