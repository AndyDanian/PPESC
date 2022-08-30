from libfock import *

class fock():
    def __init__(self, eom: list = None):
        """
        Fock object

        This object is to build or any action with fock matrix

        Args:
        ----

        eom (list): Molecular orbitals energies
        """

        if not eom or None in eom:
            print("***Warning \n\n\
                The fock object didn't recive any or nothing of the molecular orbitals energies.\n\n\
                    if you'll want calculate, you can use:\n\
                    * calculate_hf_moe: Use the wave funtion or kinetic, e->-<nucleu interactoin,\n\
                    and electron repulsion integrals")

        self._eom = eom

    ################################################################################################
    # METHODS
    ################################################################################################

    def run_hf_fock_calculate(self, intk: list = None, inten: dict = None, intee: list = None,
                            relativity_correction: bool = False, intdw: list = None, intmv: list = None,
                            mocoef: list = None, nprim: int = None, natoms: int = None, ne: int = None,
                            charge: list = None, coord: list = None, verbose: int = 0):
        """
        Run calculation Hartree--Fock molecular orbital energies

        Args:
        ----

        intk (list): 2d array with atomic kinetic integrals
        inten (dict): dictionary of 2d arrays with atomic electron--nucleus interactions integrals
        intk (list): 2d array with atomic electron repulsion integrals
        relativity_correction (bool): Activate of Massvelo, Darwin, and Spin-Orbit corrections
        intdw (list): 2d array with atomic Darwin integrals
        intmv (list): 2d array with atomic Massvelo integrals
        mocoef (list): 2d array with molecular obital coefficients
        nprim (int): primitive number
        natoms (int): atoms number
        ne (int): electrons number
        charge (list): atomic charges
        coord (list): 2d array with atomic coordinates
        verbose (int): print level

        Return:
        ------

        eom (list): 1d array with molecular orbitals energies
        """
        print_title(name = "Calculating: Fock matrix at Hartree--Fock level")

        time_start_eom = time()

        if verbose > 30:
            print_triangle_matrix(integral = mocoef, name = "Molecular Orbital Coefficients", matriz_sym = "square")

        time_start_dm = time()
        density_matrix: list = [[0.0 for zero in range(nprim)]for zero in range(nprim)]
        if ne%2 ==0:
            ne2 = int(ne/2)
        else:
            ne2 = int(ne/2) + 1
        for i in range(nprim):
            for j in range(nprim):

                for k in range(ne2):
                    density_matrix[i][j] += 2.0*mocoef[i][k]*mocoef[j][k]
        if verbose > 30:
            print_triangle_matrix(integral = density_matrix, name= "Density Matrix", matriz_sym = "sym")

        #Core Hamiltonian
        time_start_ch = time()
        hcore: list = [[0 for zero in range(nprim)]for zero in range(nprim)]
        for i in range(nprim):
            for j in range(nprim):
                ven = 0.0
                for atom_en in inten.values():
                    ven += atom_en[i][j]
                hcore[i][j] = intk[i][j] + ven

        if verbose > 30:
            print_triangle_matrix(integral = hcore, name = "Core Hamiltonian Matrix", matriz_sym = "sym")

        #Matriz G
        time_start_g = time()
        g: list = [[0 for zero in range(nprim)]for zero in range(nprim)]
        for i in range(nprim):
            for j in range(nprim):
                g[i][j] = 0.0
                for k in range(nprim):
                    for l in range(nprim):
                        g[i][j] += density_matrix[k][l]*(intee[i][j][k][l]-0.5*intee[i][l][k][j])
        if verbose > 30:
            print_triangle_matrix(integral = g, name = "G Matrix", matriz_sym = "sym")

        #Matriz Fock
        time_start_fock_ao = time()
        fock: list = [[0 for zero in range(nprim)]for zero in range(nprim)]
        if relativity_correction:
            fock_rc: list = [[0 for zero in range(nprim)]for zero in range(nprim)]

        for i  in range(nprim):
            for j  in range(nprim):
                fock[i][j] = hcore[i][j] + g[i][j]
                if relativity_correction:
                    fock_rc[i][j] = fock[i][j] + intmv[i][j] + intdw[i][j]
        if verbose > 30:
            print_triangle_matrix(integral = fock, name = "Fock Matrix in AO", matriz_sym = "square")
            if relativity_correction:
                print_triangle_matrix(integral = fock_rc,
                                name = "Fock Matrix in AO with Darwin and Massvelo Correction",
                                matriz_sym = "square")

        #FOCK
        #AO TO MO
        fock_mo: np.array = np.matmul(np.array(mocoef).T,np.matmul(np.array(fock),np.array(mocoef)))
        eom: list = [value for irow, row in enumerate(fock_mo)
                    for icol, value in enumerate(row) if irow == icol]

        if relativity_correction:
            fock_mo_rc: np.array = np.matmul(np.array(mocoef).T,np.matmul(np.array(fock_rc),np.array(mocoef)))
            eom_rc: list = [value for irow, row in enumerate(fock_mo_rc)
                    for icol, value in enumerate(row) if irow == icol]
        #Nuleu Repulsion
        time_start_te = time()
        vnn = 0.0
        distance_coordinate = [0 for i in range(3)]
        for i in range(natoms-1):
            for k in range(i + 1, natoms):
                for l in range(3):
                    distance_coordinate[l] = coord[k][l] - coord[i][l]
                    distance_coordinate[l] *= distance_coordinate[l]

                distance_magnitud = np.sqrt(sum(distance_coordinate))
                vnn += charge[i]*charge[k]/distance_magnitud

        #Electronic energy
        electronic_energy: float = 0.0
        electronic_energy_rc: float = 0.0
        for i in range(nprim):
            for j in range(nprim):
                electronic_energy += 0.5*density_matrix[i][j]*(hcore[i][j] + fock[i][j])
                if relativity_correction:
                    electronic_energy_rc += 0.5*density_matrix[i][j]*(hcore[i][j] + intdw[i][j] + intmv[i][j] + fock_rc[i][j])

        if verbose <= 10 or not verbose:
            print(f"\n Print the first 20 Hartree--Fock molecular orbitals energies: \n")
            if float(nprim/5) >= 5:
                rows: int = 5
            else:
                rows: int = int(nprim/5)
        else:
            print(f"Print all Hartree--Fock molecular orbitals energies: \n")
            rows: int = int(nprim/5)

        if nprim % 5 != 0:
            rows += 1

        for row in range(rows):
            if (row+1)*5 < nprim:
                columns = (row + 1)*5
            else:
                columns = nprim

            print(
                *[str("{:.6f}".format(eom[i])).center(14)
                for i in range(row*5, columns)],
                #end="",
            )

        if relativity_correction:
            print(f"\n\nHartree--Fock molecular orbital energies with relativities corrections: \n")
            if nprim % 5 != 0:
                rows += 1

            for row in range(rows):
                if (row+1)*5 < nprim:
                    columns = (row + 1)*5
                else:
                    columns = nprim

                print(
                    *[str("{:.4f}".format(eom_rc[i])).center(14)
                    for i in range(row*5, columns)],
                    #end="",
                )


        print("\n")
        print(40*"=")
        gap = eom[ne2]-eom[ne2-1]
        print(f"  E(LUMO):  {eom[ne2]} au")
        print(f"- E(HOMO):  {eom[ne2-1]}   au")
        print(f"----------------------------")
        print(f"    gap  :  {gap}        au")
        print(40*"=")
        print("\n")
        print(40*"=")
        print(f"\nElectronic energy: {electronic_energy}")
        print(f"Nuclear energy: {vnn}")
        print(f"Total energy (HF): {electronic_energy + vnn} \n")
        print(40*"=")
        if relativity_correction:
            print()
            print(40*"=")
            avdw: np.array = np.matmul(np.array(mocoef).T,np.matmul(np.array(intdw),np.array(mocoef)))
            avmv: np.array = np.matmul(np.array(mocoef).T,np.matmul(np.array(intmv),np.array(mocoef)))
            totaldw: float = 0.0
            totalmv: float = 0.0
            for i in range(ne2):
                totaldw += 2.0*avdw[i][i]
                totalmv += 2.0*avmv[i][i]

            total: float = totaldw + totalmv
            print(f"Darwin Correction : {totaldw:.10f} au")
            print(f"Massvelo Correction : {totalmv:.10f} au")
            print()
            print(f"Total Relativistic Corrections : {total:.10f} au ({total/(electronic_energy + vnn) * 100:.4f}%)")
            # Sumar Mv y Dw desde la construcción de la matrix de Fock, arroja el mismo resultado que
            # cuando solo sumo los valores medios de Mv y Dw al valor de la energía electrónica
            print(f"Non-Relativistic + Relativistic Corrections : {electronic_energy_rc + vnn:.10f} au")


        if verbose > 10:
            print_time(name = f"Density Matrix", delta_time = (time() - time_start_dm), tailer = False)
            print_time(name = f"Core Hamiltonian", delta_time = (time() - time_start_ch), header = False, tailer = False)
            print_time(name = f"G Matrix", delta_time = (time() - time_start_g), header = False, tailer = False)
            print_time(name = f"Fock Matrix", delta_time = (time() - time_start_fock_ao), header = False, tailer = False)
            print_time(name = f"Total Energy", delta_time = (time() - time_start_te), header = False, tailer = False)
            print_time(name = f"Molecular Orbitals Energies", delta_time = (time()-time_start_eom), header = False)

        print_title(name = "End Calculating: Fock matrix at Hartree--Fock level")
        return eom

    def calculate_hf_moe(self, wf: dict = None, intk: list = None, inten: dict = None, intee: list = None,
                        mocoef: list = None, nprim: int = None, natoms: int = None, ne: int = None,
                        charge: list = None, coord: list = None, dalton_normalization: bool = False,
                        relativity_correction: bool = False,
                        verbose: int = 0, verbose_integrals: int = 0):
        """
        Driver to calculate Hartree--Fock molecular orbital energies

        Args:
        ----

        wf (dict): dictionary with information about wave function
        intk (list): 2d array with atomic kinetic integrals
        inten (dict): dictionary of 2d arrays with atomic electron--nucleus interactions integrals
        intk (list): 2d array with atomic electron repulsion integrals
        mocoef (list): 2d array with molecular obital coefficients
        nprim (int): primitive number
        natoms (int): atoms number
        ne (int): electrons number
        charge (list): atomic charges
        coord (list): 2d array with atomic coordinates
        dalton_normalization (bool): Use dalton normalization to d, f, ...
        relativity_correction (bool): Activate of Massvelo, Darwin, and Spin-Orbit corrections
        verbose (int): print level
        verbose_integrals (int): print level for integral calculation

        Return:
        ------

        eom (list): 1d array with molecular orbitals energies

        NOTE: It's neccesary to pass list arrays to numpy arrays for time
        https://towardsdatascience.com/python-lists-are-sometimes-much-faster-than-numpy-heres-a-proof-4b3dad4653ad
        """

        start = time()

        if not wf and (not intk  or not inten or intee or not mocoef):
            raise ValueError("***Error\n\n\
                It's neccesary the wave function or kinetic, e->-<nucleu interactoin, electron \n\
                repulsion integrals, and molecular orbitals coefficient to calculate molecular\n\
                orbitals energies.")

        if wf and (not intk  or not inten or intee or not mocoef):
            calculate_integrals = eint(wf)
            print("\n\n*** Calculating: kinetic, nucpot and electron repulsion atomic integrals")

            if relativity_correction:
                print("\n\n Relativity correction Massvelo and Darwin will be calculate \n\n")
                integral_list: list = ["kinetic", "nucpot", "darwin", "massvelo"]
            else:
                integral_list: list = ["kinetic", "nucpot"]

            integrals_onebody, symmetries = calculate_integrals.integration_onebody(
                integrals_names = integral_list,
                integrals_properties = None, verbose = verbose_integrals,
                dalton_normalization = dalton_normalization)

            integrals_twobody: dict = calculate_integrals.integration_twobody(
                integrals_names = ["e2pot"], verbose = verbose_integrals,
                dalton_normalization = False
            )
            intee = integrals_twobody["e2pot"]

            # Coefficient from fock are in vector form
            mocoef: list = [list(value) for value in zip(*wf.mo_coefficients)]

            if not nprim:
                if wf._cartessian_primitive:
                    nprim: int = wf.primitives_number
                else:
                    nprim: int = wf.primitives_number_sph

            intk: list = integrals_onebody["kinetic"]
            inten: dict = {}
            for name, values in integrals_onebody.items():
                if "pot" in name:
                    inten[name] = values

            if relativity_correction:
                intdw: list = integrals_onebody["darwin"]
                intmv: list = integrals_onebody["massvelo"]

            if not charge:
                charge = wf.atomic_numbers
            if not ne:
                ne = sum(charge)
            if not coord:
                coord = wf.coordinates
            if not natoms:
                natoms = len(coord)
            #Verification mocoef would be square

        if not nprim or not ne or not natoms or not charge or not coord:
            raise ValueError("***Error\n\n\
                Information insufficient. It is necesary: primitives or atoms or\n\
                    electrons number or charges or coordinates")

        eom: list = self.run_hf_fock_calculate(intk = intk, inten = inten, intee = intee,
                                                relativity_correction = relativity_correction,
                                                intdw = intdw, intmv = intmv,
                                                mocoef = mocoef, nprim = nprim, natoms = natoms,
                                                ne = ne, charge = charge, coord = coord,
                                                verbose = verbose)

        if verbose > 10:
            print_time(name = f"Hartree-Fock Energy", delta_time = (time() - start))

        return eom

if __name__ == "__main__":
    wfn = wave_function("../tests/molden_file/HI_v2z.molden")

    print("\n Calculate MO energies used wave function \n")
    eom_values = fock()
    eom_values.calculate_hf_moe(wfn, verbose=11, verbose_integrals=1, relativity_correction=True)

# H2 STO-1G
#@    Final HF energy:              -0.160779200015
#@    Nuclear repulsion:             0.742994646761
#@    Electronic energy:            -0.903773846776

# H2 HF/STO-2G
#@    Final HF energy:              -1.095728776299
#@    Nuclear repulsion:             0.742994646761
#@    Electronic energy:            -1.838723423060
