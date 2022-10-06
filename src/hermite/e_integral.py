from typing import Union

# Current folder
from h1i import *
from h2i import *

from cto_gto_h1 import *

# Modules into sub-folder
from libint import *


class eint:
    def __init__(self, wf: wave_function):
        """
        Manages the electronic integrals calculations

        Args:
            wf (dict): dictionary with information about wave function
        """

        if not wf:
            raise ValueError(
                "*** Error \n\n There isn't information in the  wave function object."
            )

        self._wf: wave_function = wf
        # linealize arrays to calculate the integrals

        self._charge: list[float] = wf.charges
        self._coord: list[list[float]] = wf.coordinates
        self._n: int = wf.primitives_number
        self._exp: list[float] = wf.exponents
        self._center: list[int] = wf.primitives_centers
        self._lx: list[int] = wf.mlx
        self._ly: list[int] = wf.mly
        self._lz: list[int] = wf.mlz
        self._angular_moments: list[str] = wf.angular_momentums
        self._cartessian: bool = wf.cto

    ##################################################################
    # Property
    ##################################################################
    @property
    def list_1b_integrals(self) -> list[str]:
        return self._list_1b_integrals_calculated

    ##################################################################
    # METHODS
    ##################################################################
    def integration_onebody(
        self,
        integrals_names: list[str] = [],
        integrals_properties: dict[str, dict[str, list[Union[int, float]]]] = {
            "": {
                "r_gauge": [],
                "r_dipole": [],
                "spatial_symmetries": [],
                "magnetic_components": [],
                "atoms": [],
            }
        },
        gauge: list[float] = [],
        dipole: list[float] = [],
        dalton_normalization: bool = False,
        verbose: int = 0,
    ) -> None:
        """
        Driver integral one--body calculations. The integrals are
        wrote in the binary file into scratch folder in the AO2BINT.H5

        Args:
        ----
            integrals_names (list[str]): Integrals names
            integrals_properties (dict): Integrals configuations, for example, where is the
                                         gauge origen, ...
            gauge (list): Gauge for all integrals (Default)
            dipole (list): Dipole for all integrals (Default)
            Verbose: Print lever
                        0    : minimum
                        > 10 : time details
                        > 20 : atomic integrals
        """

        # Manager to write in the sratch and output file
        io = self._wf._driver_scratch
        # Manager time proces
        driver_time = self._wf._driver_time
        start: float = time()

        io.write_output(information="HERMITE: ONE BODY", type=1)

        ## Write intgral information into output ##################################
        str_integrals: str = "Integrals: "
        io.write_output(str_integrals)
        for int_name in integrals_names:
            io.write_output(
                " " * len(str_integrals) + "*" + large_name[int_name.split(" ")[0]]
            )
            if (
                int_name in integrals_properties.keys()
                and len(integrals_properties[int_name]["r_gauge"]) == 3
            ):  # in integrals_properties[int_name.split(" ")[0]].keys()):
                io.write_output(
                    " " * len(str_integrals)
                    + "   Gauge: "
                    + "{:.4f}".format(integrals_properties[int_name]["r_gauge"][0])
                    + "{:.4f}".format(integrals_properties[int_name]["r_gauge"][1])
                    + "{:.4f}".format(integrals_properties[int_name]["r_gauge"][2])
                )
        io.write_output("\n")

        if len(gauge) == 3:
            io.write_output(
                "General Gauge: "
                + "{:.4f}".format(gauge[0])
                + "{:.4f}".format(gauge[1])
                + "{:.4f}".format(gauge[2])
            )
        else:
            temp_gauge: str = (
                "{:.4f}".format(0.0) + "{:.4f}".format(0.0) + "{:.4f}".format(0.0)
            )
            io.write_output("General Gauge: " + temp_gauge)
        # END write ##############################################################

        # Verification: Is integrals implemented?
        if not integrals_names:
            raise ValueError("***Error \n\n what integral do you want?")
        else:
            for integral_name in integrals_names:
                split_integral_name = integral_name.lower().split()
                error = False
                if (
                    len(integral_name.split(" ")) == 1
                    and integral_name.lower() not in integral_symmetry.keys()
                ):
                    error = True
                elif split_integral_name[0] not in integral_symmetry.keys():
                    error = True
                if error:
                    raise ValueError(
                        f"*** Error \n\n\
                    Integral {integral_name} is not implement or the name is mistake\n\n\
                    Integrals implemented: \n\
                        {integral_symmetry.keys()}"
                    )
        # END Verification #########################################################

        # Number of atoms #########################################################
        number_atoms: int = self._wf.atom_number
        if number_atoms > 50:
            io.write_output(
                f"*** WARNING\n\n\
            System has a lot atoms ({number_atoms}), then calculate can take very much time"
            )
        # END #####################################################################

        # List varaible to save the name of the integrales calculated
        self._list_1b_integrals_calculated: list = []
        #############################################################################

        # Activation: SpinOrbit is the sum of the PSO integrals ######################
        spinorbit_integrals: bool = False
        if "spinorbit" in [
            name for int_name in integrals_names for name in int_name.split()
        ]:
            spinorbit_integrals = True
            temp_names: list = []
            activate_all_pso: bool = False
            spino_x: bool = False
            spino_y: bool = False
            spino_z: bool = False
            for name in integrals_names:
                if "spinorbit" not in name:
                    temp_names.append(name)
                else:
                    if len(name.split()) > 1:

                        if name.split()[1] == "1" or name.split()[1] == "x":
                            spino_x = True
                            self._list_1b_integrals_calculated.append("spinorbit x")
                        if name.split()[1] == "2" or name.split()[1] == "y":
                            spino_y = True
                            self._list_1b_integrals_calculated.append("spinorbit y")
                        if name.split()[1] == "3" or name.split()[1] == "z":
                            spino_z = True
                            self._list_1b_integrals_calculated.append("spinorbit z")

                        if (name.split()[1]).isalpha():
                            spatial_axe: int = [
                                i
                                for i, r in magnetic_axes.items()
                                if r == name.split()[1]
                            ][0]
                        else:
                            spatial_axe = int(name.split()[1])

                        if "pso" not in [int_name for int_name in integrals_names]:
                            temp_names += [
                                "pso " + str(spatial_axe + 1 + i * 3)
                                for i in range(number_atoms)
                                if "pso " + str(spatial_axe + 1 + i * 3)
                                not in integrals_names
                            ]

                    elif "pso" not in [int_name for int_name in integrals_names]:
                        spino_x = spino_y = spino_z = True
                        self._list_1b_integrals_calculated += [
                            "spinorbit x",
                            "spinorbit y",
                            "spinorbit z",
                        ]
                        temp_names.append("pso")
                        activate_all_pso = True
                    elif "pso" in [int_name for int_name in integrals_names]:
                        spino_x = spino_y = spino_z = True
                        self._list_1b_integrals_calculated += [
                            "spinorbit x",
                            "spinorbit y",
                            "spinorbit z",
                        ]
                        activate_all_pso = True
            if activate_all_pso:
                integrals_names = [name for name in temp_names if "pso " not in name]
            else:
                integrals_names = temp_names
        # END Spin-Orbit Integrals Activation #############################################

        # Activation: SOFIEL is the sum of the NSTCGO integrals ###########################
        sofiel_integrals: bool = False
        if "sofiel" in [
            name for int_name in integrals_names for name in int_name.split()
        ]:
            sofiel_integrals = True
            temp_names = []
            activate_all_nstcgo: bool = False
            sofiel_xx: bool = False
            sofiel_xy: bool = False
            sofiel_xz: bool = False
            sofiel_yy: bool = False
            sofiel_yz: bool = False
            sofiel_zz: bool = False
            sofiel_yx: bool = False
            sofiel_zx: bool = False
            sofiel_zy: bool = False
            for name in integrals_names:
                if "sofiel" not in name:
                    temp_names.append(name)
                else:
                    if len(name.split(" ")) > 1:
                        if name.split()[1] == "xx":
                            sofiel_xx = True
                            self._list_1b_integrals_calculated.append("sofiel xx")
                            sym_comp: int = 1
                            spatial: str = " x"
                        if name.split()[1] == "yy":
                            sofiel_yy = True
                            self._list_1b_integrals_calculated.append("sofiel yy")
                            sym_comp = 2
                            spatial = " y"
                        if name.split()[1] == "zz":
                            sofiel_zz = True
                            self._list_1b_integrals_calculated.append("sofiel zz")
                            sym_comp = 3
                            spatial = " z"
                        if name.split()[1] == "xy":
                            sofiel_xy = True
                            self._list_1b_integrals_calculated.append("sofiel xy")
                            sym_comp = 1
                            spatial = " y"
                        if name.split()[1] == "xz":
                            sofiel_xz = True
                            self._list_1b_integrals_calculated.append("sofiel xz")
                            sym_comp = 1
                            spatial = " z"
                        if name.split()[1] == "yz":
                            sofiel_yz = True
                            self._list_1b_integrals_calculated.append("sofiel yz")
                            sym_comp = 2
                            spatial = " z"
                        if name.split()[1] == "yx":
                            sofiel_yx = True
                            self._list_1b_integrals_calculated.append("sofiel yx")
                            sym_comp = 2
                            spatial = " x"
                        if name.split()[1] == "zx":
                            sofiel_zx = True
                            self._list_1b_integrals_calculated.append("sofiel zx")
                            sym_comp = 3
                            spatial = " x"
                        if name.split()[1] == "zy":
                            sofiel_zy = True
                            self._list_1b_integrals_calculated.append("sofiel zy")
                            sym_comp = 3
                            spatial = " y"

                        if "nstcgo" not in [int_name for int_name in integrals_names]:
                            temp_names += [
                                "nstcgo " + str(sym_comp + i * 3) + spatial
                                for i in range(number_atoms)
                                if "nstcgo " + str(sym_comp + i * 3) + spatial
                                not in integrals_names
                            ]
                    elif "nstcgo" not in [int_name for int_name in integrals_names]:
                        sofiel_xx = (
                            sofiel_yy
                        ) = (
                            sofiel_zz
                        ) = (
                            sofiel_xy
                        ) = (
                            sofiel_xz
                        ) = sofiel_yz = sofiel_yx = sofiel_zx = sofiel_zy = True
                        self._list_1b_integrals_calculated += [
                            "sofiel xx",
                            "sofiel yy",
                            "sofiel zz",
                            "sofiel xy",
                            "sofiel xz",
                            "sofiel yz",
                            "sofiel yx",
                            "sofiel zx",
                            "sofiel zy",
                        ]
                        temp_names.append("nstcgo")
                        activate_all_nstcgo = True
                    elif "nstcgo" in [int_name for int_name in integrals_names]:
                        sofiel_xx = (
                            sofiel_yy
                        ) = (
                            sofiel_zz
                        ) = (
                            sofiel_xy
                        ) = (
                            sofiel_xz
                        ) = sofiel_yz = sofiel_yx = sofiel_zx = sofiel_zy = True
                        self._list_1b_integrals_calculated += [
                            "sofiel xx",
                            "sofiel yy",
                            "sofiel zz",
                            "sofiel xy",
                            "sofiel xz",
                            "sofiel yz",
                            "sofiel yx",
                            "sofiel zx",
                            "sofiel zy",
                        ]

                        activate_all_nstcgo = True
            if activate_all_nstcgo:
                integrals_names = [name for name in temp_names if "nstcgo " not in name]
            else:
                integrals_names = temp_names
        # END SOFIEL integrals activation ####################################################

        # Remove repet integrals #############################################################
        temp_integrals_names: list = []
        for count, int_name in enumerate(integrals_names):
            if int_name.split()[0] in integrals_names[count + 1 :]:
                continue
            if count > 0 and int_name.split()[0] in temp_integrals_names:
                continue
            temp_integrals_names.append(int_name)
        integrals_names = temp_integrals_names
        # END Remove repet integrals ##########################################################

        # * Start 1B Integral Calculation ******************************************************
        symmetries: dict[str, str] = {}
        for int_name in integrals_names:
            integrals: dict[str, list[float]] = {}

            if len(int_name.split(" ")) > 1:
                integral_name = int_name.lower().split()[0]
            else:
                integral_name = int_name

            # Check definition of integrals properties into integrals_properties
            if integral_name in integrals_properties.keys():
                if integrals_properties[integral_name]["r_gauge"] == 3:
                    r_gauge = integrals_properties[integral_name]["r_gauge"]
                if integrals_properties[integral_name]["r_dipole"] == 3:
                    r_dipole = integrals_properties[integral_name]["r_dipole"]
                if len(integrals_properties[integral_name]["magnetic_components"]) > 0:
                    magnetic_components = [
                        x - 1
                        for x in integrals_properties[integral_name][
                            "magnetic_components"
                        ]
                    ]
                if len(integrals_properties[integral_name]["spatial_symmetries"]) > 0:
                    spatial_symmetries = [
                        s - 1
                        for s in integrals_properties[integral_name][
                            "spatial_symmetries"
                        ]
                    ]
                if len(integrals_properties[integral_name]["atoms"]) > 0:
                    atoms = [
                        a - 1 for a in integrals_properties[integral_name]["atoms"]
                    ]
            else:
                (
                    r_gauge,
                    r_dipole,
                    magnetic_components,
                    spatial_symmetries,
                    atoms,
                ) = integral_1b_parameters(
                    atoms_number=number_atoms,
                    integral_name=int_name,
                    gauge=gauge,
                    dipole=dipole,
                )
            # END Check integrals properties ########################################################

            if (
                spatial_symmetry[integral_name.lower()] == 0
                and magnetic[integral_name.lower()] == 0
            ):

                if integral_name.lower() in [
                    "overlap",
                    "darwin",
                    "massvelo",
                    "kinetic",
                ]:

                    # Build Name == integral name #############################
                    integral_label: str = integral_name.lower()
                    # END Build Name ##########################################

                    # Verification: Integrals is calculated, already ###########
                    # - Save integral name
                    self._list_1b_integrals_calculated.append(integral_label)
                    if io._hermite_ao1b_binary.exists() and io.binary(
                        file=io._hermite_ao1b_binary, label=integral_label, io="f"
                    ):
                        continue
                    # END Verfication ##########################################

                    symmetries[integral_label] = integral_symmetry[
                        integral_name.lower()
                    ]
                    integrals[integral_label] = h1i(
                        # Default
                        charge=self._charge,
                        coord=self._coord,
                        exp=self._exp,
                        center=self._center,
                        lx=self._lx,
                        ly=self._ly,
                        lz=self._lz,
                        name=integral_name,
                        verbose=verbose,
                        dalton_normalization=dalton_normalization,
                        driver_time=driver_time,
                    )
                elif integral_name.lower() in ["nucpot", "fc"]:

                    for atom in atoms:
                        # Build Name according atom index ################
                        integral_label = integral_name.lower() + " " + str(atom + 1)
                        # End Build Name #################################

                        # Verification: Integrals is calculated, already ###########
                        # - Save integral name
                        self._list_1b_integrals_calculated.append(integral_label)
                        if io._hermite_ao1b_binary.exists() and io.binary(
                            file=io._hermite_ao1b_binary, label=integral_label, io="f"
                        ):
                            continue
                        # END Verification ##########################################

                        symmetries[integral_label] = integral_symmetry[
                            integral_name.lower()
                        ]
                        integrals[integral_label] = h1i(
                            # Default
                            charge=self._charge,
                            coord=self._coord,
                            exp=self._exp,
                            center=self._center,
                            lx=self._lx,
                            ly=self._ly,
                            lz=self._lz,
                            name=integral_name,
                            verbose=verbose,
                            dalton_normalization=dalton_normalization,
                            driver_time=driver_time,
                            # Special information
                            atom=atom,
                        )

            elif (
                spatial_symmetry[integral_name.lower()] == 0
                and magnetic[integral_name.lower()] == 1
            ):

                for b_i in magnetic_components:

                    # Build name of the integral according magnetic or spatial component ###
                    if type(b_i) == int:
                        if integral_name.lower() == "laplacian":
                            integral_label = (
                                integral_name.lower() + " " + spatial_components[b_i]
                            )
                        else:
                            integral_label = (
                                integral_name.lower() + " " + magnetic_axes[b_i]
                            )
                        magnetic_xyz: int = b_i
                    else:
                        integral_label = integral_name.lower() + " " + str(b_i)
                        if integral_name.lower() == "laplacian":
                            magnetic_xyz = list(spatial_components.keys())[
                                list(spatial_components.values()).index(b_i)
                            ]
                        else:
                            magnetic_xyz = list(magnetic_axes.keys())[
                                list(magnetic_axes.values()).index(b_i)
                            ]
                    # END Build Name ####################################################

                    # Verification: Integrals is calculated, already #####################
                    # - Save integral name
                    self._list_1b_integrals_calculated.append(integral_label)
                    if io._hermite_ao1b_binary.exists() and io.binary(
                        file=io._hermite_ao1b_binary, label=integral_label, io="f"
                    ):
                        continue
                    # END Verification ###################################################

                    symmetries[integral_label] = integral_symmetry[
                        integral_name.lower()
                    ]

                    integrals[integral_label] = h1i(
                        # Default
                        charge=self._charge,
                        coord=self._coord,
                        exp=self._exp,
                        center=self._center,
                        lx=self._lx,
                        ly=self._ly,
                        lz=self._lz,
                        name=integral_name,
                        verbose=verbose,
                        dalton_normalization=dalton_normalization,
                        driver_time=driver_time,
                        # Special information
                        magnetic_xyz=magnetic_xyz,
                        r_gauge=r_gauge,
                        r_dipole=r_dipole,
                    )

            elif (
                spatial_symmetry[integral_name.lower()] == 1
                and magnetic[integral_name.lower()] == 0
            ):

                for spatial_i in spatial_symmetries:

                    # Selection of coordinate x, y, z for spatial symmetry #######
                    coordinate: int = int(spatial_i) - 3 * int(spatial_i / 3)
                    atom_index: int = int(spatial_i / 3)
                    # END Selection ################################################

                    # Verification atom index ######################################
                    if atom_index >= number_atoms:
                        raise ValueError(
                            f"***Error \n\n\
                            atom {atom_index} doesn't exist"
                        )
                    # END Verification ############################################

                    # Spatial component ############################################
                    if coordinate >= 0 and coordinate <= 2:
                        spatial_component: int = coordinate
                    else:
                        raise ValueError(
                            f"*** Error\n\n \
                            spatial component doesn't exist, {spatial_component}"
                        )
                    # END Spatial Component ########################################

                    # Build Name ##################################################
                    integral_label = str(
                        integral_name.lower() + " " + str(spatial_i + 1)
                    )
                    # END Build Name ###############################################

                    # Verification: Integrals is calculated, already #################
                    # - Save integral name
                    self._list_1b_integrals_calculated.append(integral_label)
                    if io._hermite_ao1b_binary.exists() and io.binary(
                        file=io._hermite_ao1b_binary, label=integral_label, io="f"
                    ):
                        continue
                    # END Verification ###############################################

                    symmetries[integral_label] = integral_symmetry[
                        integral_name.lower()
                    ]
                    integrals[integral_label] = h1i(
                        # Default
                        charge=self._charge,
                        coord=self._coord,
                        exp=self._exp,
                        center=self._center,
                        lx=self._lx,
                        ly=self._ly,
                        lz=self._lz,
                        name=integral_name,
                        verbose=verbose,
                        dalton_normalization=dalton_normalization,
                        driver_time=driver_time,
                        # Special information
                        spatial_sym=spatial_component,
                        atom=atom_index,
                    )

            elif (
                spatial_symmetry[integral_name.lower()] == 1
                and magnetic[integral_name.lower()] == 1
            ):

                #! Faltaría implementar cuando solo se da la componente spatial or magnetic, para calcular
                #! esta con todas las componentes de la otra

                for spatial_i in spatial_symmetries:

                    # Selection of coordinate x, y, z for spatial symmetry#########
                    coordinate = int(spatial_i) - 3 * int(spatial_i / 3)
                    atom_index = int(spatial_i / 3)
                    # END Selection Coordinate ####################################

                    # Verification: atom index ####################################
                    if atom_index >= number_atoms:
                        raise ValueError(
                            f"***Error \n\n\
                            atom {atom_index} doesn't exist"
                        )
                    # END Verification ############################################

                    # Spatial Component ###########################################
                    if coordinate >= 0 and coordinate <= 2:
                        spatial_component = coordinate
                    else:
                        raise ValueError(
                            "*** Error\n\n \
                            spatial sym doesn't exist"
                        )
                    # End Spatial Component ########################################

                    for b_i in magnetic_components:

                        if type(spatial_i) == int:
                            integral_label = (
                                integral_name.lower()
                                + " "
                                + str(spatial_i + 1)
                                + " "
                                + magnetic_axes[b_i]
                            )
                            magnetic_xyz = int(b_i)
                        else:
                            integral_label = (
                                integral_name.lower()
                                + " "
                                + str(spatial_i + 1)
                                + " "
                                + str(b_i)
                            )
                            magnetic_xyz = list(magnetic_axes.keys())[
                                list(magnetic_axes.values()).index(b_i)
                            ]

                        # Veridication: Integrals is calculated, already #####################
                        # - Save integral name
                        self._list_1b_integrals_calculated.append(integral_label)
                        if io._hermite_ao1b_binary.exists() and io.binary(
                            file=io._hermite_ao1b_binary, label=integral_label, io="f"
                        ):
                            continue
                        # END Verification #####################################################

                        symmetries[integral_label] = integral_symmetry[
                            integral_name.lower()
                        ]
                        integrals[integral_label] = h1i(
                            # Default
                            charge=self._charge,
                            coord=self._coord,
                            exp=self._exp,
                            center=self._center,
                            lx=self._lx,
                            ly=self._ly,
                            lz=self._lz,
                            name=integral_name,
                            verbose=verbose,
                            dalton_normalization=dalton_normalization,
                            driver_time=driver_time,
                            # Special information
                            magnetic_xyz=magnetic_xyz,
                            spatial_sym=spatial_component,
                            atom=atom_index,
                            r_gauge=r_gauge,
                        )
            # * END 1B Integral Calculation ****************************************************************

            # Transform integrals from cto to sph
            integrals_matrix = {}
            if not self._cartessian:
                start_cto: float = time()
                for label, integral in integrals.items():
                    integrals_matrix[label] = cto_gto_h1(
                        vector_to_matrix(self._n, integral, symmetries[label]),
                        np.array(self._angular_moments),
                    )
                time_cto: float = time() - start_cto
            else:
                for label, integral in integrals.items():
                    integrals_matrix[label] = vector_to_matrix(
                        self._n, integral, symmetries[label]
                    )
            ## Write in binary file
            for label, integral in integrals_matrix.items():
                io.binary(
                    file=io._hermite_ao1b_binary, dictionary={label: integral}, io="a"
                )

        # * Calculation: SpinOrbit Calculation and Write integrals in AO1BINT **************
        if spinorbit_integrals:
            spin_orbit(
                integrals=io,
                driver_time=driver_time,
                list_spino=[
                    spino_x,
                    spino_y,
                    spino_z,
                ],
                charge=self._charge,
                number_atoms=number_atoms,
                nprim=self._wf.primitives_number,
                verbose=verbose,
            )
        # *END SPINORBITCalculation ******************************************************
        # * Calculation: SOFIEL Calculation and Write integrals in AO1BINT ****************
        if sofiel_integrals:
            sofiel(
                integrals=io,
                number_atoms=number_atoms,
                charge=self._charge,
                nprim=self._wf.primitives_number,
                list_sofiel=[
                    sofiel_xx,
                    sofiel_xy,
                    sofiel_xz,
                    sofiel_yx,
                    sofiel_yy,
                    sofiel_yz,
                    sofiel_zx,
                    sofiel_zy,
                    sofiel_zz,
                ],
                driver_time=driver_time,
                verbose=verbose,
            )
        # *END SOFIEL Calculation *********************************************************
        ## Write in output the size of AO1BINT.H5 in bytes
        io.write_output(
            information=io._hermite_ao1b_binary.name,
            type=3,
            size_file=io._hermite_ao1b_binary.stat().st_size,
        )
        # Write integrals in the output file
        if verbose > 20:
            io.write_output(type=9)
        # Time
        if verbose > 10:
            driver_time.add_name_delta_time(
                name=f"One--Body CTOs--GTOs", delta_time=time_cto
            )

        driver_time.add_name_delta_time(
            name="Hermite Calculation", delta_time=(time() - start)
        )

        # Write time into output file
        io.write_output("\n")
        io.write_output(type=2, drv_time=driver_time)
        driver_time.reset

        io.write_output(information="END HERMITE: ONE BODY", type=1)

    def integration_twobody(
        self,
        integrals_names: list = [],
        verbose: int = 0,
        dalton_normalization: bool = False,
    ) -> None:
        """
        Driver to calculation the two--body atomic integrals. The integrals are
        wrote in the binary file into scratch folder in the AO2BINT.H5.

        Implemented:
            repulsion integrals
        """

        ## Instace external objects
        # - Scratch
        io = self._wf._driver_scratch
        io.write_output(information="HERMITE: TWO BODY", type=1)
        # - Diver Time
        driver_time = self._wf._driver_time
        ##
        start = time()

        if len(integrals_names) == 0:
            integral_name: str = "e2pot"

        integrals_2_cart: dict[str, np.ndarray] = {}
        integral_name = integrals_names[0]

        # * Start 1B Integral Calculation ******************************************************
        integrals_2_cart[integral_name.lower()] = h2i(
            # Default
            coord=self._coord,
            exp=self._exp,
            center=self._center,
            lx=self._lx,
            ly=self._ly,
            lz=self._lz,
            name=integral_name,
            verbose=verbose,
            dalton_normalization=dalton_normalization,
            driver_time=driver_time,
        )
        # * END 1B Integral Calculation ********************************************************

        integrals_two_body: dict = {}
        # Cartessian to Spherical
        if not self._cartessian:
            start_cto: float = time()
            for label, integral in integrals_2_cart.items():
                integrals_two_body[label] = cart2sph(
                    intcart=np.asfortranarray(integral),
                    angular=np.array(self._angular_moments),
                    nprim=len(self._angular_moments),
                    ncar=self._wf.primitives_number_car,
                    nsph=self._wf.primitives_number_sph,
                )
            driver_time.add_name_delta_time(
                name=f"Two--Body CTOs--GTOs", delta_time=time() - start_cto
            )
        else:
            for label, integral in integrals_2_cart.items():
                integrals_two_body[label] = np.array(integral)

        ## Write in binary file two body atomic integrals
        io.binary(file=io._hermite_ao2b_binary, dictionary=integrals_two_body, io="a")
        ## Write in output the AO2BINT.H5 size in bytes
        io.write_output(
            information=io._hermite_ao2b_binary.name,
            type=3,
            size_file=io._hermite_ao2b_binary.stat().st_size,
        )

        if verbose > 100:
            if not self._cartessian:
                information: str = "Two--body integrals with gto--primitives"
            else:
                information = "Two--body integrals with cto--primitives"
            io.write_output(information=information, type=1, title_type=1)
            # Write integrals in the output file
            start_print: float = time()
            io.write_output(type=10)
            driver_time.add_name_delta_time(
                name="Print two integrals", delta_time=time() - start_print
            )

        driver_time.add_name_delta_time(
            name="Hermite Calculation", delta_time=(time() - start)
        )

        # Write time into output file
        io.write_output("\n")
        io.write_output(type=2, drv_time=driver_time)
        driver_time.reset

        io.write_output(information="END HERMITE: TWO BODY", type=1)


if __name__ == "__main__":
    wf = wave_function(
        "../tests/molden_file/LiH.molden",
        scratch_path="/home1/scratch",
        job_folder="160922134451",
    )
    s = eint(wf)
    one = True
    if one:
        s.integration_onebody(
            integrals_names=["nucpot", "darwin", "fc", "spinorbit", "sofiel"],
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
            verbose=31,
            dalton_normalization=False,
        )
    else:
        integrals = s.integration_twobody(
            ["e2pot"], verbose=11, dalton_normalization=False
        )
