# Current folder
from h1i import *
from h2i import *
# Modules into sub-folder
from libint import *
from cto_gto_h1 import *
from cto_gto_h2 import *


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
            self._angular_moments,
            self._cartessian,
        ) = linealize_array_wf(wf)

    ##################################################################
    # METHODS
    ##################################################################

    def integration_onebody(
        self, integrals_names: list = None, integrals_properties: dict = None, output: int = 0,
        dalton_normalization: bool = False
    ):

        if not integrals_names:
            raise ValueError("***Error \n\n what integral do you want?")
        else:
            for integral_name in integrals_names:
                if integral_name.lower() not in integral_symmetry.keys():
                    raise ValueError(f"*** Error \n\n\
                    Integral {integral_name} is not implement or the name is mistake\n\n\
                    Integrals implemented: \n\
                        {integral_symmetry.keys()}")

        # output dictionaries
        integrals: dict = {}
        symmetries: dict = {}

        # Default values to integrals properties
        r_gauge: list = [0.0,0.0,0.0]
        r_dipole: list = [0.0,0.0,0.0]
        atoms: list = [atom for atom in range(len(self._coord))]
        magnetic_components: list = [0,1,2]
        spatial_symmetries: list = [symmetry for symmetry in range(np.size(np.array(self._coord)))]
        number_atoms: int =  len(self._coord[:][0])

        if len(atoms) > 50:
            print(f"*** WARNING\n\n\
            System has a lot atoms ({len(atoms)}), then calculate can take very much time")

        for integral_name in integrals_names:

        # Check definition of integrals properties into integrals_properties
            if integrals_properties:
                if integral_name in integrals_properties.keys():
                    if "r_gauge" in integrals_properties[integral_name].keys():
                        r_gauge = integrals_properties[integral_name]["r_gauge"]
                    if "r_dipole" in integrals_properties[integral_name].keys():
                        r_dipole = integrals_properties[integral_name]["r_dipole"]
                    if "magnetic_components" in integrals_properties[integral_name].keys():
                        magnetic_components = integrals_properties[integral_name]["magnetic_components"]
                    if "spatial_symmetries" in integrals_properties[integral_name].keys():
                        spatial_symmetries = integrals_properties[integral_name]["spatial_symmetries"]
                    if "atoms" in integrals_properties[integral_name].keys():
                        atoms = integrals_properties[integral_name]["atoms"]

            if spatial_symmetry[integral_name.lower()] == 0 and magnetic[integral_name.lower()] == 0:

                if integral_name.lower() in ["overlap", "darwin", "kinetic"]:

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
                        output = output,
                        dalton_normalization = dalton_normalization
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
                            output = output,
                            dalton_normalization = dalton_normalization,
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
                        output = output,
                        dalton_normalization = dalton_normalization,
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
                        output = output,
                        dalton_normalization = dalton_normalization,
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
                            output = output,
                            dalton_normalization = dalton_normalization,
                            # Special information
                            magnetic_xyz = magnetic_xyz,
                            spatial_sym = spatial_component,
                            atom = atom,
                            r_gauge = r_gauge
                        )

        # Print integral
        integrals_matrix = {}
        if not self._cartessian:
            for label, integral in integrals.items():
                integrals_matrix[label] = cto_gto_h1(np.array(vector_to_matrix(len(self._exp), integral, symmetries[label])),
                        np.array(self._angular_moments))
            if output > 20:
                print("="*80,"\n One--body integrals with gto--primitives\n","="*80)
                print_matriz_integrated(integrals = integrals_matrix, symmetries = symmetries)
        else:
            for label, integral in integrals.items():
                integrals_matrix[label] = np.array(vector_to_matrix(len(self._exp), integral, symmetries[label]))
            if output > 20:
                print("="*80,"\n One--body integrals with cto--primitives\n",80*"=")
                print_matriz_integrated(integrals = integrals_matrix, symmetries = symmetries)


        return integrals_matrix, symmetries


    def integration_twobody(
        self, integrals_names: list = None,output: int = 0,
        dalton_normalization: bool = False
    ):
        """
        Driver to calculation the two--body atomic integrals

        Implemented:
            repulsion integrals
        """

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
            output = output,
            dalton_normalization = dalton_normalization
        )

        integrals_two_body = {}
        if output > 100 and self._cartessian:
            for label, integral in integrals_2_cart.items():
                integrals_two_body[label] = np.array(integral)
            print("="*80,"\nTwo--body integrals with cto--primitives\n","="*80)
            print(integrals_two_body["e2pot"])

        if not self._cartessian:
            for label, integral in integrals_2_cart.items():
                integrals_two_body[label] = cto_gto_h2(np.array(integral),
                        np.array(self._angular_moments))

        if output > 100:
            print("="*80,"\nTwo--body integrals with gto--primitives\n","="*80)
            print(integrals_two_body["e2pot"])

        return integrals_two_body



if __name__ == "__main__":
    from libint import *
    wfn = wave_function("../io/LiH_sd.molden")

    s = eint(wfn.build_wfn_array())

    one = False
    if one:
        integrals, symmetries = s.integration_onebody(["kinetic"],
                    {
                    "nucpot":{"atoms":[0]},
                    "angmom":{"magnetic_components":[0, 1, 2], "r_gauge":[0.0, 0.0, 1.404552358700]},
                    "sd":{"spatial_symmetries":[0,1,2,3,4,5], "magnetic_components":[0,1,2]},
                    "fc":{"atoms":[0,1]},
                    "nelfld":{"spatial_symmetries":[0,1,2,3,4,5]},
                    "diplen":{"r_dipole":[0.0,0.0,0.0],"magnetic_components":[0,1,2]},
                    "dipvel":{"magnetic_components":[0,1,2]},
                    "pso":{"spatial_symmetries":[0,1,2,3,4,5]},
                    "nstcgo":{"spatial_symmetries":[0,1,2,3,4,5],"magnetic_components":[0,1,2], "r_gauge":[0.0, 0.0, 1.404552358700]},
                    "dnske":{"spatial_symmetries":[0,1,2,3,4,5],"magnetic_components":[0,1,2], "r_gauge":[0.0, 0.0, 1.404552358700]},
                    "psoke":{"spatial_symmetries":[0,1,2,3,4,5]},
                    "psooz":{"spatial_symmetries":[0,1,2,3,4,5],"magnetic_components":[0,1,2], "r_gauge":[0.0, 0.0, 1.404552358700]},
                    "ozke":{"magnetic_components":[0,1,2], "r_gauge":[0.0, 0.0, 1.404552358700]},
                    },
                    21, dalton_normalization=False)

    integrals = s.integration_twobody(["e2pot"], output=11, dalton_normalization=False)
