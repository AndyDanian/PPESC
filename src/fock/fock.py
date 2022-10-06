from libfock import *


class fock:
    def __init__(self, wf: wave_function):
        """
        Fock object

        This object is to build or any action with fock matrix

        Args:
        ----
        wf (dict): dictionary with information about wave function
        eom (list): Molecular orbitals energies
        """

        self._wf: wave_function = wf

    ################################################################################################
    # Properties
    ################################################################################################
    @property
    def hf_moe(self) -> list:
        return self.hf_eom_calculated

    ################################################################################################
    # METHODS
    ################################################################################################

    def run_hf_fock_calculate(
        self,
        io: scratch,
        driver_time: drv_time,
        intk: np.ndarray,
        inten: dict[str, np.ndarray],
        intee: np.ndarray,
        intdw: np.ndarray,
        intmv: np.ndarray,
        mocoef: list,
        charge: list,
        coord: list,
        nprim: int,
        natoms: int,
        ne: int,
        relativity_correction: bool = False,
        verbose: int = 0,
    ):
        """
        Run calculation Hartree--Fock molecular orbital energies

        Args:
        ----
        io (object:scratch): Driver to driver the output and binary files
        time_object (drv_time): Manage time calculation
        intk (np.ndarray): 2d array with atomic kinetic integrals
        inten (dict[str, np.ndarray]): dictionary of 2d arrays with atomic electron--nucleus interactions integrals
        intk (np.ndarray): 2d array with atomic electron repulsion integrals
        intdw (np.ndarray): 2d array with atomic Darwin integrals
        intmv (np.ndarray): 2d array with atomic Massvelo integrals
        mocoef (list): 2d array with molecular obital coefficients
        nprim (int): primitive number
        natoms (int): atoms number
        ne (int): electrons number
        charge (list): atomic charges
        coord (list): 2d array with atomic coordinates
        relativity_correction (bool): Activate of Massvelo, Darwin, and Spin-Orbit corrections
        verbose (int): print level

        Return:
        ------

        eom (list): 1d array with molecular orbitals energies
        """
        io.write_output(
            information="Calculating: Fock matrix at Hartree--Fock level", type=1
        )

        time_start_eom: float = time()

        if verbose > 30:
            io.write_output(
                type=9,
                direct=True,
                dictionary={"Molecular Orbital Coefficients": mocoef},
            )

        time_start_dm: float = time()
        density_matrix: np.ndarray = np.zeros((nprim, nprim), dtype=float)

        if ne % 2 == 0:
            ne2 = int(ne / 2)
        else:
            ne2 = int(ne / 2) + 1

        for i in range(nprim):
            for j in range(nprim):
                for k in range(ne2):
                    density_matrix[i, j] += 2.0 * mocoef[i][k] * mocoef[j][k]
        if verbose > 30:
            io.write_output(
                type=9, direct=True, dictionary={"Density Matrix": density_matrix}
            )

        # Core Hamiltonian
        time_start_ch: float = time()
        hcore: np.ndarray = np.zeros((nprim, nprim), dtype=float)
        for i in range(nprim):
            for j in range(nprim):
                ven: float = 0.0
                for atom_en in inten.values():
                    ven += atom_en[i, j]
                hcore[i, j] = intk[i, j] + ven
        if verbose > 30:
            io.write_output(
                type=9, direct=True, dictionary={"Core Hamiltonian Matrix": hcore}
            )

        # Matriz G
        time_start_g: float = time()
        g: np.ndarray = np.zeros((nprim, nprim), dtype=float)
        for i in range(nprim):
            for j in range(nprim):
                g[i, j] = 0.0
                for k in range(nprim):
                    for l in range(nprim):
                        g[i, j] += density_matrix[k, l] * (
                            intee[i][j][k][l] - 0.5 * intee[i][l][k][j]
                        )
        if verbose > 30:
            io.write_output(type=9, direct=True, dictionary={"G Matrix": g})

        # Matriz Fock
        time_start_fock_ao: float = time()
        fock: np.ndarray = np.zeros((nprim, nprim), dtype=float)
        for i in range(nprim):
            for j in range(nprim):
                fock[i, j] = hcore[i, j] + g[i, j]
        if verbose > 30:
            io.write_output(type=9, direct=True, dictionary={"Fock Matrix in AO": fock})

        # FOCK
        # AO TO MO
        fock_mo: np.ndarray = np.matmul(
            np.array(mocoef).T, np.matmul(np.array(fock), np.array(mocoef))
        )
        eom: list = [fock_mo[i][i] for i in range(nprim)]
        self.hf_eom_calculated: list = eom

        # Nuleu Repulsion
        time_start_te: float = time()
        vnn: float = 0.0
        distance_coordinate: list[float] = [0.0 for i in range(3)]
        for i in range(natoms - 1):
            for k in range(i + 1, natoms):
                for l in range(3):
                    distance_coordinate[l] = coord[k][l] - coord[i][l]
                    distance_coordinate[l] *= distance_coordinate[l]

                distance_magnitud = np.sqrt(sum(distance_coordinate))
                vnn += charge[i] * charge[k] / distance_magnitud

        # Electronic energy
        electronic_energy: float = 0.0
        for i in range(nprim):
            for j in range(nprim):
                electronic_energy += (
                    0.5 * density_matrix[i][j] * (hcore[i][j] + fock[i][j])
                )

        if verbose <= 10 or not verbose:
            io.write_output(
                f"\n Print the first 20 Hartree--Fock molecular orbitals energies: \n"
            )
            if float(nprim / 5) >= 5:
                rows: int = 5
            else:
                rows = int(nprim / 5)
        else:
            io.write_output(f"Print all Hartree--Fock molecular orbitals energies: \n")
            rows = int(nprim / 5)

        if nprim % 5 != 0:
            rows += 1

        for row in range(rows):
            if (row + 1) * 5 < nprim:
                columns: int = (row + 1) * 5
            else:
                columns = nprim

            row_strings: str = " "
            for i in range(row * 5, columns):

                if abs(eom[i]) > 9999.0:
                    formate: str = "{:.6e}"
                else:
                    formate = "{:.6f}"

                row_strings += formate.format(eom[i]).center(14) + " "
            io.write_output(row_strings)

        io.write_output("\n")
        gap: float = eom[ne2] - eom[ne2 - 1]
        io.write_output(f"  E(LUMO): " + f"{eom[ne2]:4f}".center(24) + " au")
        io.write_output(f"- E(HOMO): " + f"{eom[ne2-1]:4f}".center(24) + " au")
        io.write_output(f"-" * 40)
        io.write_output("    gap  : " + f"{gap:4f}".center(24) + " au")

        io.write_output("\n")
        io.write_output(40 * "=")
        io.write_output(f"\nElectronic energy: {electronic_energy}")
        io.write_output(f"Nuclear energy: {vnn}")
        io.write_output(f"Total energy (HF): {electronic_energy + vnn} \n")
        io.write_output(40 * "=")
        # Relativity correction for the energy
        if relativity_correction:

            intmv_mo: np.ndarray = np.matmul(
                np.array(mocoef).T, np.matmul(np.array(intmv), np.array(mocoef))
            )
            intdw_mo: np.ndarray = np.matmul(
                np.array(mocoef).T, np.matmul(np.array(intdw), np.array(mocoef))
            )
            if verbose <= 10 or not verbose:
                io.write_output(
                    f"\nPrint the first 20 Hartree--Fock molecular orbitals energies with\n\
                                    relativities corrections: \n"
                )
                if float(nprim / 3) >= 3:
                    rows = 3
                else:
                    rows = int(nprim / 3)
            else:
                io.write_output(
                    f"\n\nHartree--Fock molecular orbitals energies with relativities corrections: \n"
                )
                rows = int(nprim / 3)

            if nprim % 3 != 0:
                rows += 1

            for row in range(rows):
                if (row + 1) * 3 < nprim:
                    columns = (row + 1) * 3
                else:
                    columns = nprim

                row_strings = " "

                for i in range(row * 3, columns):
                    energy_relativity: float = eom[i] + intmv_mo[i][i] + intdw_mo[i][i]
                    percent: float = (
                        (intmv_mo[i][i] + intdw_mo[i][i])
                        / abs(eom[i] + intmv_mo[i][i] + intdw_mo[i][i])
                        * 100
                    )
                    if abs(energy_relativity) > 9999.0:
                        formate = "{:.6e}({:.3f}%)"
                    else:
                        formate = "{:.6f}({:.3f}%)"

                    row_strings += (
                        formate.format(energy_relativity, percent).center(28) + " "
                    )
                io.write_output(row_strings)

            io.write_output(" ")
            io.write_output(
                "Note: Between parethesis is the relativity correction in porcentage\n\
                            (formule: relativity correction/ABS(Total)*100)"
            )
            io.write_output(" ")

            for i in range(nprim):
                eom[i] += intmv_mo[i][i] + intdw_mo[i][i]
            self.hf_eom_calculated = eom

            gap = eom[ne2] - eom[ne2 - 1]
            io.write_output(f"  E(LUMO): " + f"{eom[ne2]:4f}".center(24) + " au")
            io.write_output(f"- E(HOMO): " + f"{eom[ne2-1]:4f}".center(24) + " au")
            io.write_output(f"-" * 40)
            io.write_output("    gap  : " + f"{gap:4f}".center(24) + " au")

            count: int = 0
            for i in range(ne2, nprim):
                if eom[i] < 0.0:
                    count += 1
            if count > 0:
                io.write_output("*" * 70)
                io.write_output("***WARNING***")
                if count == 1:
                    io.write_output(f"The last virtual orbital energy is negative.\n\n")
                    io.write_output(
                        f"Then, there is a exponent in the basis set very high, in general,"
                    )
                else:
                    io.write_output(
                        f"The last {count} virtual orbitals energies are negative.\n\n"
                    )
                    io.write_output(
                        f"Then, there are exponents in the basis set very high, in general,"
                    )
                io.write_output(f"they are associated with S atomic orbitals.")
                io.write_output("*" * 70)

            io.write_output(" ")
            io.write_output(80 * "=")
            totaldw: float = 0.0
            totalmv: float = 0.0
            for i in range(ne2):
                totaldw += 2.0 * intdw_mo[i][i]
                totalmv += 2.0 * intmv_mo[i][i]

            total: float = totaldw + totalmv
            io.write_output(f"Darwin Correction : {totaldw:.10f} au")
            io.write_output(f"Massvelo Correction : {totalmv:.10f} au")
            io.write_output(" ")
            io.write_output(
                f"Total Relativistic Corrections : {total:.10f} au ({total/(electronic_energy + vnn) * 100:.4f}%)"
            )
            # Sumar Mv y Dw desde la construcción de la matrix de Fock, arroja el mismo resultado que
            # cuando solo sumo los valores medios de Mv y Dw al valor de la energía electrónica
            io.write_output(
                f"Non-Relativistic + Relativistic Corrections : {electronic_energy + vnn + total:.10f} au"
            )
            io.write_output(80 * "=")

        if verbose > 10:
            driver_time.add_name_delta_time(
                name=f"Density Matrix", delta_time=(time() - time_start_dm)
            )
            driver_time.add_name_delta_time(
                name=f"Core Hamiltonian", delta_time=(time() - time_start_ch)
            )
            driver_time.add_name_delta_time(
                name=f"G Matrix", delta_time=(time() - time_start_g)
            )
            driver_time.add_name_delta_time(
                name=f"Fock Matrix", delta_time=(time() - time_start_fock_ao)
            )
            driver_time.add_name_delta_time(
                name=f"Total Energy", delta_time=(time() - time_start_te)
            )
            driver_time.add_name_delta_time(
                name=f"Molecular Orbitals Energies",
                delta_time=(time() - time_start_eom),
            )

        io.write_output(
            information="End Calculating: Fock matrix at Hartree--Fock level", type=1
        )

        return eom

    def calculate_hf_moe(
        self,
        relativity_correction: bool = False,
        dalton_normalization: bool = False,
        verbose: int = 0,
        verbose_integrals: int = 0,
    ) -> list[float]:
        """
        Driver to calculate Hartree--Fock molecular orbital energies

        Args:
        ----
        dalton_normalization (bool): Use dalton normalization to d, f, ...
        relativity_correction (bool): Activate of Massvelo, Darwin, and Spin-Orbit corrections
        verbose (int): print level
        verbose_integrals (int): print level for integral calculation
            > 10: Time more details
            > 30: Print arrays

        Return:
        ------

        eom (list): 1d array with molecular orbitals energies

        NOTE: It's neccesary to pass list arrays to numpy arrays for time
        https://towardsdatascience.com/python-lists-are-sometimes-much-faster-than-numpy-heres-a-proof-4b3dad4653ad
        """

        ## Instance external objects
        # - Scratch
        io = self._wf._driver_scratch
        # - Driver Time
        driver_time = self._wf._driver_time
        ##
        start: float = time()

        if not self._wf:
            raise ValueError(
                "***Error\n\n\
                It's neccesary the wave function to calculate molecular\n\
                orbitals energies."
            )

        calculate_integrals = eint(self._wf)
        io.write_output(
            "\n\n*** Calculating: kinetic, nucpot and electron repulsion atomic integrals"
        )
        if relativity_correction:
            io.write_output(
                "    Relativity correction Massvelo and Darwin will be calculate \n\n"
            )
            integral_list: list = ["kinetic", "nucpot", "darwin", "massvelo"]
        else:
            integral_list = ["kinetic", "nucpot"]

        calculate_integrals.integration_onebody(
            integrals_names=integral_list,
            verbose=verbose_integrals,
            dalton_normalization=dalton_normalization,
        )

        calculate_integrals.integration_twobody(
            integrals_names=["e2pot"],
            verbose=verbose_integrals,
            dalton_normalization=False,
        )
        intee: np.ndarray = io.binary(
            file=io._hermite_ao2b_binary, io="r", label="e2pot"
        )
        # Coefficient from fock are in vector form
        mocoef = [list(value) for value in zip(*self._wf.mo_coefficients)]

        if self._wf._cartessian_primitive:
            nprim: int = self._wf.primitives_number
        else:
            nprim = self._wf.primitives_number_sph

        intk: np.ndarray = io.binary(
            file=io._hermite_ao1b_binary, io="r", label="kinetic"
        )
        inten: dict[str, np.ndarray] = {}
        for a in range(self._wf.atom_number):
            name = f"pot {a + 1}"
            inten[name] = io.binary(
                file=io._hermite_ao1b_binary, io="r", label="kinetic"
            )

        if relativity_correction:
            intdw: np.ndarray = io.binary(
                file=io._hermite_ao1b_binary, io="r", label="darwin"
            )
            intmv: np.ndarray = io.binary(
                file=io._hermite_ao1b_binary, io="r", label="massvelo"
            )

        # Atomic charge
        charge: list = self._wf.atomic_numbers
        # Electrons number
        ne: int = sum(charge)
        # Atomic coordinates
        coord: list[list[float]] = self._wf.coordinates
        # Atom number
        natoms: int = self._wf.atom_number

        if not nprim or not ne or not natoms or not charge or not coord:
            raise ValueError(
                "***Error\n\n\
                Information insufficient. It is necesary: primitives or atoms or\n\
                    electrons number or charges or coordinates"
            )

        eom: list = self.run_hf_fock_calculate(
            io=io,
            driver_time=driver_time,
            intk=intk,
            inten=inten,
            intee=intee,
            relativity_correction=relativity_correction,
            intdw=intdw,
            intmv=intmv,
            mocoef=mocoef,
            nprim=nprim,
            charge=charge,
            coord=coord,
            natoms=natoms,
            ne=ne,
            verbose=verbose,
        )

        driver_time.add_name_delta_time(
            name=f"Hartree-Fock Energy", delta_time=(time() - start)
        )
        io.write_output(type=2, drv_time=driver_time)
        driver_time.reset

        return eom


if __name__ == "__main__":
    wfn = wave_function(
        "../tests/molden_file/H2.molden",
        scratch_path="/home1/scratch",
        job_folder="160922134451",
    )

    print("\n Calculate MO energies used wave function \n")
    eom_values = fock(wfn)
    eom_values.calculate_hf_moe(
        verbose=11, verbose_integrals=1, relativity_correction=True
    )

# H2 STO-1G
# @    Final HF energy:              -0.160779200015
# @    Nuclear repulsion:             0.742994646761
# @    Electronic energy:            -0.903773846776

# H2 HF/STO-2G
# @    Final HF energy:              -1.095728776299
# @    Nuclear repulsion:             0.742994646761
# @    Electronic energy:            -1.838723423060
