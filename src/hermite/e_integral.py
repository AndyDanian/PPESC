# Current folder
from h1i import *
from h2i import *
from cto_gto_h1 import *
from cto_gto_h2 import *
# Modules into sub-folder
from libint import *


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

        self._wf = wf
        # linealize arrays to calculate the integrals

        self._charge = wf.charges
        self._coord  = wf.coordinates
        self._n      = wf.primitives_number
        self._exp    = wf.exponents
        self._center = wf.primitives_centers
        self._lx     = wf.mlx
        self._ly     = wf.mly
        self._lz     = wf.mlz
        self._angular_moments = wf.angular_momentums
        self._cartessian = wf.cto

    ##################################################################
    # METHODS
    ##################################################################

    def integration_onebody(
        self, integrals_names: list = None, integrals_properties: dict = None, verbose: int = 0,
        gauge: list  = None, dipole: list = None, dalton_normalization: bool = False
    ):

        if verbose >= 0:
            print_title(name = "HERMITE: ONE BODY")

        if verbose > 10:
            driver_time = drv_time()
            start = time()
        else:
            driver_time = None

        if not integrals_names:
            raise ValueError("***Error \n\n what integral do you want?")
        else:
            for integral_name in integrals_names:
                split_integral_name = integral_name.lower().split()
                error = False
                if len(integral_name.split(" ")) == 1 and integral_name.lower() not in integral_symmetry.keys():
                    error = True
                elif split_integral_name[0] not in integral_symmetry.keys():
                    error = True
                if error:
                    raise ValueError(f"*** Error \n\n\
                    Integral {integral_name} is not implement or the name is mistake\n\n\
                    Integrals implemented: \n\
                        {integral_symmetry.keys()}")

        # verbose dictionaries
        integrals: dict = {}
        symmetries: dict = {}

        number_atoms: int =  self._wf.atom_number
        if number_atoms > 50:
            print(f"*** WARNING\n\n\
            System has a lot atoms ({len(atoms)}), then calculate can take very much time")

        if gauge is not None:
            r_gauge = gauge
        if dipole is not None:
            r_dipole = dipole

        ## SpinOrbit is the sum of the PSO integrals
        spinorbit_integrals: bool = False
        if "spinorbit" in [name for int_name in integrals_names for name in int_name.split()]:
            spinorbit_integrals: bool = True
            temp_names: list = []
            activate_all_pso: bool = False
            for name in integrals_names:
                if "spinorbit" not in name:
                    temp_names.append(name)
                else:
                    if len(name.split(" ")) > 1 and "pso" not in [int_name for int_name in integrals_names]:
                        if name.split()[1] == "1": spino_x: bool = True
                        if name.split()[1] == "2": spino_y: bool = True
                        if name.split()[1] == "3": spino_z: bool = True
                        temp_names += ["pso " + str(int(name.lower().split()[1]) + i*3)
                                        for i in range(number_atoms)
                                        if "pso " + str(int(name.lower().split()[1]) + i*3) not in integrals_names]
                    elif "pso" not in [int_name for int_name in integrals_names]:
                        spino_x = spino_y = spino_z = True
                        temp_names.append("pso")
                        activate_all_pso = True
            if activate_all_pso:
                integrals_names = [name for name in temp_names if "pso " not in name]

        for int_name in integrals_names:
            if len(int_name.split(" ")) > 1:
                integral_name = int_name.lower().split()[0]
            else:
                integral_name = int_name

        # Check definition of integrals properties into integrals_properties
            if integrals_properties:
                if integral_name in integrals_properties.keys():
                    if "r_gauge" in integrals_properties[integral_name].keys():
                        r_gauge = integrals_properties[integral_name]["r_gauge"]
                    if "r_dipole" in integrals_properties[integral_name].keys():
                        r_dipole = integrals_properties[integral_name]["r_dipole"]
                    if "magnetic_components" in integrals_properties[integral_name].keys():
                        magnetic_components = [x - 1 for x in integrals_properties[integral_name]["magnetic_components"]]
                    if "spatial_symmetries" in integrals_properties[integral_name].keys():
                        spatial_symmetries = [s - 1 for s  in integrals_properties[integral_name]["spatial_symmetries"]]
                    if "atoms" in integrals_properties[integral_name].keys():
                        atoms = [a - 1 for a in integrals_properties[integral_name]["atoms"]]
            else: # When is indicated in the name the magnetic or symmetry name or atom
                r_gauge, r_dipole, magnetic_components, spatial_symmetries, atoms =\
                    integral_1b_parameters(atoms_number = number_atoms, integral_name = int_name)

            if spatial_symmetry[integral_name.lower()] == 0 and magnetic[integral_name.lower()] == 0:

                if integral_name.lower() in ["overlap", "darwin", "massvelo", "kinetic"]:

                    symmetries[integral_name.lower()] = integral_symmetry[integral_name.lower()]
                    integrals[integral_name.lower()] = h1i(
                        #Default
                        charge = self._charge,
                        coord = self._coord,
                        exp = self._exp,
                        center = self._center,
                        lx = self._lx,
                        ly = self._ly,
                        lz = self._lz,
                        name = integral_name,
                        verbose = verbose,
                        dalton_normalization = dalton_normalization,
                        driver_time = driver_time
                    )
                elif integral_name.lower() in ["nucpot", "fc"]:

                    for atom in atoms:
                        symmetries[
                            integral_name.lower() + " " + str(atom + 1)
                            ] = integral_symmetry[integral_name.lower()]
                        integrals[
                            integral_name.lower() + " " + str(atom + 1)
                            ] = h1i(
                            # Default
                            charge = self._charge,
                            coord = self._coord,
                            exp = self._exp,
                            center = self._center,
                            lx = self._lx,
                            ly = self._ly,
                            lz = self._lz,
                            name = integral_name,
                            verbose = verbose,
                            dalton_normalization = dalton_normalization,
                            driver_time = driver_time,
                            # Special information
                            atom = atom,
                        )

            elif spatial_symmetry[integral_name.lower()] == 0 and magnetic[integral_name.lower()] == 1:

                for b_i in magnetic_components:

                    if type(b_i) == int:
                        integral_label: str = (integral_name.lower() + " " + magnetic_axes[b_i])
                        magnetic_xyz: int = b_i
                    else:
                        integral_label: str = (
                            integral_name.lower() + " " + b_i)
                        magnetic_xyz: int = (list(magnetic_axes.keys())
                                            [list(magnetic_axes.values()).index(b_i)])

                    symmetries[integral_label] = integral_symmetry[integral_name.lower()]

                    integrals[integral_label] = h1i(
                        #Default
                        charge = self._charge,
                        coord = self._coord,
                        exp = self._exp,
                        center = self._center,
                        lx = self._lx,
                        ly = self._ly,
                        lz = self._lz,
                        name = integral_name,
                        verbose = verbose,
                        dalton_normalization = dalton_normalization,
                        driver_time = driver_time,
                        # Special information
                        magnetic_xyz = magnetic_xyz,
                        r_gauge = r_gauge,
                        r_dipole = r_dipole
                    )

            elif spatial_symmetry[integral_name.lower()] == 1 and magnetic[integral_name.lower()] == 0:

                for spatial_i in spatial_symmetries:

                    # Selection of coordinate x, y, z for spatial symmetry
                    coordinate: int = spatial_i - 3 * int(spatial_i/3)
                    atom: int = int(spatial_i/3)

                    if atom >= number_atoms:
                        raise ValueError(f"***Error \n\n\
                            atom {atom} doesn't exist")

                    if coordinate == 0:
                        spatial_component: int = 0
                    elif coordinate == 1:
                        spatial_component: int = 1
                    elif coordinate == 2:
                        spatial_component: int = 2
                    else:
                        raise ValueError("*** Error\n\n \
                            spatial component doesn't exist, {spatial_component}")

                    integral_label: str = str(
                        integral_name.lower() + " " + str(spatial_i + 1))

                    symmetries[integral_label] = integral_symmetry[integral_name.lower()]
                    integrals[integral_label] = h1i(
                        # Default
                        charge = self._charge,
                        coord = self._coord,
                        exp = self._exp,
                        center = self._center,
                        lx = self._lx,
                        ly = self._ly,
                        lz = self._lz,
                        name = integral_name,
                        verbose = verbose,
                        dalton_normalization = dalton_normalization,
                        driver_time = driver_time,
                        # Special information
                        spatial_sym = spatial_component,
                        atom = atom,
                    )

            elif spatial_symmetry[integral_name.lower()] == 1 and magnetic[integral_name.lower()] == 1:

                for spatial_i in spatial_symmetries:

                    # Selection of coordinate x, y, z for spatial symmetry
                    coordinate: int = spatial_i - 3 * int(spatial_i/3)
                    atom: int = int(spatial_i/3)

                    if atom >= number_atoms:
                        raise ValueError(f"***Error \n\n\
                            atom {atom} doesn't exist")

                    if coordinate == 0:
                        spatial_component: int = 0
                    elif coordinate == 1:
                        spatial_component: int = 1
                    elif coordinate == 2:
                        spatial_component: int = 2
                    else:
                        raise ValueError("*** Error\n\n \
                            spatial sym doesn't exist")

                    for b_i in magnetic_components:

                        if type(spatial_i) == int:
                            integral_label: str = (integral_name.lower() + " " +
                            str(spatial_i + 1) + " " + magnetic_axes[b_i])
                            magnetic_xyz: int = b_i
                        else:
                            integral_label: str = (integral_name.lower() + " " +
                            str(spatial_i + 1) + " "  + b_i)
                            magnetic_xyz: int = (list(magnetic_axes.keys())
                            [list(magnetic_axes.values()).index(b_i)])

                        symmetries[integral_label] = integral_symmetry[integral_name.lower()]

                        integrals[integral_label] = h1i(
                            # Default
                            charge = self._charge,
                            coord = self._coord,
                            exp = self._exp,
                            center = self._center,
                            lx = self._lx,
                            ly = self._ly,
                            lz = self._lz,
                            name = integral_name,
                            verbose = verbose,
                            dalton_normalization = dalton_normalization,
                            driver_time = driver_time,
                            # Special information
                            magnetic_xyz = magnetic_xyz,
                            spatial_sym = spatial_component,
                            atom = atom,
                            r_gauge = r_gauge
                        )

        ### SpinOrbit Calculation
        if spinorbit_integrals:
            temp_spinorbit_integrals: dict = {}
            if spino_x:
                symmetries["spinorbit x"] = "antisym"
                for a in range(number_atoms):
                    if a == 0:
                        temp_spinorbit_integrals["spinorbit x"] = [
                                        self._charge[a] * value for value in integrals["pso " + str(1 + a*3)]
                                        ]
                    else:
                        temp_spinorbit_integrals["spinorbit x"] = [
                                    old + self._charge[a] * new
                                    for old, new in
                                        zip(temp_spinorbit_integrals["spinorbit x"], integrals["pso " + str(1 + a*3)])]
                driver_time.add_name_delta_time(name = "Spin-Orbit X AO", delta_time = (time() - start))
            if spino_y:
                symmetries["spinorbit y"] = "antisym"
                for a in range(number_atoms):
                    if a == 0:
                        temp_spinorbit_integrals["spinorbit y"] = [
                                        self._charge[a] * value for value in integrals["pso " + str(2 + a*3)]
                                        ]
                    else:
                        temp_spinorbit_integrals["spinorbit y"] = [
                                    old + self._charge[a] * new
                                    for old, new in
                                        zip(temp_spinorbit_integrals["spinorbit y"], integrals["pso " + str(2 + a*3)])]
                driver_time.add_name_delta_time(name = "Spin-Orbit Y AO", delta_time = (time() - start))
            if spino_z:
                symmetries["spinorbit z"] = "antisym"
                for a in range(number_atoms):
                    if a == 0:
                        temp_spinorbit_integrals["spinorbit z"] = [
                                        self._charge[a] * value for value in integrals["pso " + str(3 + a*3)]
                                        ]
                    else:
                        temp_spinorbit_integrals["spinorbit z"] = [
                                    old + self._charge[a] * new
                                    for old, new in
                                        zip(temp_spinorbit_integrals["spinorbit z"], integrals["pso " + str(3 + a*3)])]
                driver_time.add_name_delta_time(name = "Spin-Orbit Z AO", delta_time = (time() - start))
            integrals.update(temp_spinorbit_integrals)

        # Print integral
        integrals_matrix = {}
        if not self._cartessian:
            for label, integral in integrals.items():
                integrals_matrix[label] = cto_gto_h1(np.array(vector_to_matrix(self._n, integral, symmetries[label])),
                        np.array(self._angular_moments), driver_time = driver_time)
            if verbose > 20:
                print_title(name = "One--body integrals with gto--primitives")
                print_matriz_integrated(integrals = integrals_matrix, symmetries = symmetries)
        else:
            for label, integral in integrals.items():
                integrals_matrix[label] = np.array(vector_to_matrix(self._n, integral, symmetries[label]))
            if verbose > 20:
                print_title(name = "One--body integrals with cto--primitives")
                print_matriz_integrated(integrals = integrals_matrix, symmetries = symmetries)

        if verbose > 10:
            driver_time.add_name_delta_time(name = "Hermite Calculation", delta_time = (time() - start))
            driver_time.printing()

        if verbose >= 0:
            print_title(name = f"END HERMITE: ONE BODY")

        return integrals_matrix, symmetries


    def integration_twobody(
        self, integrals_names: list = None, verbose: int = 0,
        dalton_normalization: bool = False
    ):
        """
        Driver to calculation the two--body atomic integrals

        Implemented:
            repulsion integrals
        """

        if verbose >= 0:
            print_title(name = "HERMITE: TWO BODY")

        if verbose > 10:
            driver_time = drv_time()
            start = time()
        else:
            driver_time = None


        if integrals_names == None:
            integral_name: str = "e2pot"

        integrals_2_cart: dict = {}
        integral_name: str = integrals_names[0]

        integrals_2_cart[integral_name.lower()] = h2i(
            #Default
            coord = self._coord,
            exp = self._exp,
            center = self._center,
            lx = self._lx,
            ly = self._ly,
            lz = self._lz,
            name = integral_name,
            verbose = verbose,
            dalton_normalization = dalton_normalization,
            driver_time = driver_time
        )

        integrals_two_body = {}
        if verbose > 100 and self._cartessian:
            for label, integral in integrals_2_cart.items():
                integrals_two_body[label] = np.array(integral)
            print("="*80,"\nTwo--body integrals with cto--primitives\n","="*80)
            print(integrals_two_body["e2pot"])

        # Cartessian to Spherical
        if not self._cartessian:
            for label, integral in integrals_2_cart.items():
                integrals_two_body[label] = cto_gto_h2(np.array(integral),
                        np.array(self._angular_moments), driver_time = driver_time)

        if verbose > 100:
            print("="*80,"\nTwo--body integrals with gto--primitives\n","="*80)
            print(integrals_two_body["e2pot"])

        if verbose > 10:
            driver_time.add_name_delta_time(name = "Hermite Calculation", delta_time = (time() - start))
            driver_time.printing()

        if verbose >= 0:
            print_title(name = f"END HERMITE: TWO BODY")

        return integrals_two_body



if __name__ == "__main__":
    wf = wave_function("../tests/molden_file/LiH.molden")
    s = eint(wf)
    one = True
    if one:
        integrals, symmetries = s.integration_onebody(integrals_names = ["spinorbit", "pso 3"],
                    # {
                    # "nucpot":{"atoms":[0]},
                    # "angmom":{"magnetic_components":[0, 1, 2], "r_gauge":[0.0, 0.0, 1.404552358700]},
                    # "sd":{"spatial_symmetries":[0,1,2,3,4,5], "magnetic_components":[0,1,2]},
                    # "fc":{"atoms":[0,1]},
                    # "nelfld":{"spatial_symmetries":[0,1,2,3,4,5]},
                    # "diplen":{"r_dipole":[0.0,0.0,0.0],"magnetic_components":[0,1,2]},
                    # "dipvel":{"magnetic_components":[0,1,2]},
                    # "pso":{"spatial_symmetries":[0,1,2,3,4,5]},
                    # "nstcgo":{"spatial_symmetries":[0,1,2,3,4,5],"magnetic_components":[0,1,2], "r_gauge":[0.0, 0.0, 1.404552358700]},
                    # "dnske":{"spatial_symmetries":[0,1,2,3,4,5],"magnetic_components":[0,1,2], "r_gauge":[0.0, 0.0, 1.404552358700]},
                    # "psoke":{"spatial_symmetries":[0,1,2,3,4,5]},
                    # "psooz":{"spatial_symmetries":[0,1,2,3,4,5],"magnetic_components":[0,1,2], "r_gauge":[0.0, 0.0, 1.404552358700]},
                    # "ozke":{"magnetic_components":[0,1,2], "r_gauge":[0.0, 0.0, 1.404552358700]},
                    # },
                    verbose = 21, dalton_normalization=False)
    else:
        integrals = s.integration_twobody(["e2pot"], verbose=11, dalton_normalization=False)
