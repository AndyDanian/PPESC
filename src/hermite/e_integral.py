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
        gaugeo: list  = None, dipole: list = None, dalton_normalization: bool = False
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

        ## SpinOrbit is the sum of the PSO integrals
        spinorbit_integrals: bool = False
        if "spinorbit" in [name for int_name in integrals_names for name in int_name.split()]:
            spinorbit_integrals: bool = True
            temp_names: list = []
            activate_all_pso: bool = False
            spino_x = spino_y = spino_z = False
            old_integrals_names: list = integrals_names
            for name in integrals_names:
                if "spinorbit" not in name:
                    temp_names.append(name)
                else:
                    if len(name.split(" ")) > 1 and "pso" not in [int_name for int_name in integrals_names]:
                        if name.split()[1] == "1" or name.split()[1] == "x": spino_x: bool = True
                        if name.split()[1] == "2" or name.split()[1] == "y": spino_y: bool = True
                        if name.split()[1] == "3" or name.split()[1] == "z": spino_z: bool = True

                        if (name.split()[1]).isalpha():
                            spatial_axe: int = [i for i ,r in magnetic_axes.items() if r == name.split()[1]][0]
                        else:
                            spatial_axe: int = int(name.split()[1])
                        temp_names += ["pso " + str(spatial_axe + 1 + i*3)
                                        for i in range(number_atoms)
                                        if "pso " + str(spatial_axe + 1 + i*3) not in integrals_names]
                    elif "pso" not in [int_name for int_name in integrals_names]:
                        spino_x = spino_y = spino_z = True
                        temp_names.append("pso")
                        activate_all_pso = True
                    elif "pso" in [int_name for int_name in integrals_names]:
                        spino_x = spino_y = spino_z = True
                        activate_all_pso = True
            if activate_all_pso:
                integrals_names = [name for name in temp_names if "pso " not in name]
            else:
                integrals_names = temp_names


        ## SOFIEL is the sum of the NSTCGO integrals
        sofiel_integrals: bool = False
        if "sofiel" in [name for int_name in integrals_names for name in int_name.split()]:
            sofiel_integrals: bool = True
            temp_names: list = []
            activate_all_nstcgo: bool = False
            sofiel_xx = sofiel_xy = sofiel_xz = sofiel_yy = sofiel_yz = sofiel_zz =\
                sofiel_yx = sofiel_zx = sofiel_zy = False
            old_integrals_names: list = integrals_names
            for name in integrals_names:
                if "sofiel" not in name:
                    temp_names.append(name)
                else:
                    if len(name.split(" ")) > 1 and "nstcgo" not in [int_name for int_name in integrals_names]:
                        if name.split()[1] == "xx":
                            sofiel_xx: bool = True
                            sym_comp: int = 1
                            spatial: str = " x"
                        if name.split()[1] == "yy":
                            sofiel_yy: bool = True
                            sym_comp: int = 2
                            spatial: str = " y"
                        if name.split()[1] == "zz":
                            sofiel_zz: bool = True
                            sym_comp: int = 3
                            spatial: str = " z"
                        if name.split()[1] == "xy":
                            sofiel_xy: bool = True
                            sym_comp: int = 1
                            spatial: str = " y"
                        if name.split()[1] == "xz":
                            sofiel_xz: bool = True
                            sym_comp: int = 1
                            spatial: str = " z"
                        if name.split()[1] == "yz":
                            sofiel_yz: bool = True
                            sym_comp: int = 2
                            spatial: str = " z"
                        if name.split()[1] == "yx":
                            sofiel_yx: bool = True
                            sym_comp: int = 2
                            spatial: str = " x"
                        if name.split()[1] == "zx":
                            sofiel_zx: bool = True
                            sym_comp: int = 3
                            spatial: str = " x"
                        if name.split()[1] == "zy":
                            sofiel_zy: bool = True
                            sym_comp: int = 3
                            spatial: str = " y"
                        temp_names += ["nstcgo " + str(sym_comp + i*3) + spatial
                                        for i in range(number_atoms)
                                        if "nstcgo " + str(sym_comp + i*3) + spatial not in integrals_names]
                    elif "nstcgo" not in [int_name for int_name in integrals_names]:
                        sofiel_xx = sofiel_yy = sofiel_zz = sofiel_xy = sofiel_xz = sofiel_yz =\
                            sofiel_yx = sofiel_zx = sofiel_zy = True
                        temp_names.append("nstcgo")
                        activate_all_nstcgo = True
                    elif "nstcgo" in [int_name for int_name in integrals_names]:
                        sofiel_xx = sofiel_yy = sofiel_zz = sofiel_xy = sofiel_xz = sofiel_yz =\
                            sofiel_yx = sofiel_zx = sofiel_zy = True
                        activate_all_nstcgo = True
            if activate_all_nstcgo:
                integrals_names = [name for name in temp_names if "nstcgo " not in name]
            else:
                integrals_names = temp_names

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
            else:
                r_gauge, r_dipole, magnetic_components, spatial_symmetries, atoms =\
                                integral_1b_parameters(atoms_number = number_atoms, integral_name = int_name,
                                                        gaugeo = gaugeo, dipole = dipole)

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
                        if integral_name.lower() == "laplacian":
                            integral_label: str = (integral_name.lower() + " " + spatial_components[b_i])
                        else:
                            integral_label: str = (integral_name.lower() + " " + magnetic_axes[b_i])
                        magnetic_xyz: int = b_i
                    else:
                        integral_label: str = (
                            integral_name.lower() + " " + b_i)
                        if integral_name.lower() == "laplacian":
                            magnetic_xyz: int = (list(spatial_components.keys())
                                            [list(spatial_components.values()).index(b_i)])
                        else:
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
                        raise ValueError(f"*** Error\n\n \
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

                #! Faltaría implementar cuando solo se da la componente spatial or magnetic, para calcular
                #! esta con todas las componentes de la otra

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

        # Transform integrals from cto to sph
        integrals_matrix = {}
        if not self._cartessian:
            start_cto: float = time()
            for label, integral in integrals.items():
                integrals_matrix[label] = cto_gto_h1(np.array(vector_to_matrix(self._n, integral, symmetries[label])),
                        np.array(self._angular_moments))
            time_cto: float = time() - start_cto
        else:
            for label, integral in integrals.items():
                integrals_matrix[label] = np.array(vector_to_matrix(self._n, integral, symmetries[label]))

        ### SpinOrbit Calculation
        if spinorbit_integrals:
            temp_integrals_matrix: dict = {}
            so_integrals, so_symmetries = spin_orbit(integrals = integrals_matrix, number_atoms = number_atoms,
                                                    charge = self._charge, nprim = self._wf.primitives_number,
                                                    spino_x = spino_x, spino_y = spino_y, spino_z = spino_z,
                                                    driver_time = driver_time, verbose = verbose)
            # Errase pso
            temp_integrals_matrix = {name: value for name, value in integrals_matrix.items()
                                    if name.lower().split()[0] in old_integrals_names}
            integrals_matrix = so_integrals
            integrals_matrix.update(temp_integrals_matrix)
            symmetries.update(so_symmetries)
        ### SOFIEL Calculation
        if sofiel_integrals:
            temp_integrals_matrix: dict = {}
            sf_integrals,  sf_symmetries = sofiel(integrals = integrals_matrix, number_atoms = number_atoms,
                                                charge = self._charge, nprim = self._wf.primitives_number,
                                                sofiel_xx = sofiel_xx, sofiel_yy = sofiel_yy, sofiel_zz = sofiel_zz,
                                                sofiel_xy = sofiel_xy, sofiel_xz = sofiel_xz, sofiel_yz = sofiel_yz,
                                                sofiel_yx = sofiel_yx, sofiel_zx = sofiel_zx, sofiel_zy = sofiel_zy,
                                                driver_time = driver_time, verbose = verbose)
            # Errase nstcgo
            temp_integrals_matrix = {name: value for name, value in integrals_matrix.items()
                                    if name.lower().split()[0] in old_integrals_names}
            integrals_matrix = sf_integrals
            integrals_matrix.update(temp_integrals_matrix)
            symmetries.update(sf_symmetries)

        # Print Integrals
        if not self._cartessian and verbose > 20:
            print_title(name = "One--body integrals with gto--primitives")
            print_matriz_integrated(integrals = integrals_matrix, symmetries = symmetries)
        elif verbose > 20:
            print_title(name = "One--body integrals with cto--primitives")
            print_matriz_integrated(integrals = integrals_matrix, symmetries = symmetries)

        # Time
        if verbose > 10:
            driver_time.add_name_delta_time(name = f"One--Body CTOs--GTOs", delta_time = time_cto)
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
            start_cto: float = time()
            for label, integral in integrals_2_cart.items():
                integrals_two_body[label] = cto_gto_h2(np.array(integral),
                        np.array(self._angular_moments))
            time_cto = time() - start_cto
        if verbose > 100:
            print("="*80,"\nTwo--body integrals with gto--primitives\n","="*80)
            print(integrals_two_body["e2pot"])

        if verbose > 10:
            driver_time.add_name_delta_time(name = f"Two--Body CTOs--GTOs", delta_time = time_cto)
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
        integrals, symmetries = s.integration_onebody(integrals_names = ["nucpot","darwin"],
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
        integrals = s.integration_twobody(["e2pot"], verbose=10, dalton_normalization=False)
