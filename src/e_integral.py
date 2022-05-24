from lib import *
from h1i import *
from wave_function import *

nucleu: dict = {"overlap": 0, "nucpot": 1, "kinetic": 0, "angmom": 0, "sd": 0, "fc": 1, "darwin": 0,
"massvelo": 0, "nelfld": 0, "diplen": 0, "dipvel": 0, "pso": 0, "nstcgo": 0, "dnske": 0, "psoke": 0,
"psooz": 0, "ozke": 0}
esp_sym: dict = {"overlap": 0, "nucpot": 0, "kinetic": 0, "angmom": 0, "sd": 1, "fc": 0, "darwin": 0,
"massvelo": 0, "nelfld": 1, "diplen": 0, "dipvel": 0, "pso": 1, "nstcgo": 1, "dnske": 1, "psoke": 1,
"psooz": 1, "ozke": 0}
magnetic_r: dict = {"overlap": 0, "nucpot": 0, "kinetic": 0, "angmom": 1, "sd": 1, "fc": 0, "darwin": 0,
"massvelo": 0, "nelfld": 0, "diplen": 1, "dipvel": 1, "pso": 0, "nstcgo": 1, "dnske": 1, "psoke": 0,
"psooz": 1, "ozke": 1}
integral_symmetry: dict = {"overlap": "sym", "nucpot": "sym", "kinetic": "sym", "angmom": "antisym",
"sd": "sym", "fc": "sym", "darwin": "sym", "massvelo": "sym", "nelfld": "sym", "diplen": "sym", 
"dipvel": "antisym", "pso": "antisym", "nstcgo": "sym", "dnske": "sym", "psoke": "square",
"psooz": "square", "ozke": "antisym"}

magnetic_components: dict = {0:"x", 1:"y", 2:"z"}

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
        self, names: list = None, properties: list = None, output: int = None
    ):

        if not names:
            raise ValueError("***Error \n\n what integral do you want?")
        else:
            for name in names:
                if name.lower() not in integral_symmetry.keys():
                    raise ValueError(f"*** Error \n\n\
                    Integral name is not implement or the name is mistake\n\n\
                    Integrals implemented: \n\
                        {integral_symmetry.keys()}")

        integrals = {}
        symmetry = {}

        for name in names:
            if nucleu[name.lower()] == 0 and esp_sym[name.lower()] == 0 and magnetic_r[name.lower()] == 0:
                symmetry[str(name)] = integral_symmetry[name.lower()]

                integrals[str(name)] = h1i(
                    charge = self._charge,
                    coord = self._coord,
                    exp = self._exp,
                    center = self._center,
                    lx = self._lx,
                    ly = self._ly,
                    lz = self._lz,
                    name = name,
                    output = output
                )

            if nucleu[name.lower()] == 1 and esp_sym[name.lower()] == 0  and magnetic_r[name.lower()] == 0:
                for atom in properties[name]["atoms"]:
                    symmetry[
                        name.lower() + " " + str(atom + 1)
                    ] = integral_symmetry[name.lower()]

                    integrals[name.lower() + " " + str(atom + 1)] = h1i(
                        charge = self._charge,
                        coord = self._coord,
                        exp = self._exp,
                        center = self._center,
                        lx = self._lx,
                        ly = self._ly,
                        lz = self._lz,
                        name = name,
                        output = output,
                        atom = atom,
                    )

            if  nucleu[name.lower()] == 0 and magnetic_r[name.lower()] == 1 and esp_sym[name.lower()] == 0:
                
                if "rdipole" not in properties[name].keys():
                    properties[name]["rdipole"] = None

                if "gauge" not in properties[name].keys():
                    properties[name]["gauge"] = None


                for m_component in properties[name]["magnetic"]:

                    if type(m_component) == int:
                        integral_name: str = (name.lower() + " " 
                        + magnetic_components[m_component])
                        magnetic_xyz: int = m_component
                    else:
                        integral_name: str = (name.lower() + " " 
                        + m_component)
                        magnetic_xyz: int = (list(magnetic_components.keys())
                                            [list(magnetic_components.values()).index(m_component)])

                    symmetry[
                        integral_name
                    ] = integral_symmetry[name.lower()]

                    integrals[integral_name] = h1i(
                        charge = self._charge,
                        coord = self._coord,
                        exp = self._exp,
                        center = self._center,
                        lx = self._lx,
                        ly = self._ly,
                        lz = self._lz,
                        name = name,
                        output = output,
                        magnetic_xyz = magnetic_xyz,
                        gauge = properties[name]["gauge"],
                        rdipole = properties[name]["rdipole"]
                    )

            if nucleu[name.lower()] == 0 and magnetic_r[name.lower()] == 0 and esp_sym[name.lower()] == 1:
                number_atoms: int =  len(self._coord[:][0])

                for spatial in properties[name]["spatial"]:

                    # Selection of coordinate x, y, z for spatial symmetry
                    coordinate: int = spatial - 3 * int(spatial/3)
                    atom: int = int(spatial/3)
                    
                    if atom >= number_atoms:
                        raise ValueError(f"***Error \n\n\
                            atom {atom} doesn't exist") 
    
                    if coordinate == 0:
                        spatial_sym: int = 0
                    elif coordinate == 1:
                        spatial_sym: int = 1
                    elif coordinate == 2:
                        spatial_sym: int = 2
                    else:
                        raise ValueError("*** Error\n\n \
                            spatial sym doesn't exist")

                    if type(spatial) == int:
                        integral_name: str = (name.lower() + " " +
                        str(spatial + 1))
                    else:
                        integral_name: str = (name.lower() + " " +
                        str(spatial + 1))

                    symmetry[
                        integral_name
                    ] = integral_symmetry[name.lower()]

                    integrals[integral_name] = h1i(
                        charge = self._charge,
                        coord = self._coord,
                        exp = self._exp,
                        center = self._center,
                        lx = self._lx,
                        ly = self._ly,
                        lz = self._lz,
                        name = name,
                        output = output,
                        spatial_sym = spatial_sym,
                        atom = atom
                    )

            if nucleu[name.lower()] == 0 and magnetic_r[name.lower()] == 1 and esp_sym[name.lower()] == 1:
                number_atoms: int =  len(self._coord[:][0])

                if "gauge" not in properties[name].keys():
                    properties[name]["gauge"] = None

                for spatial in properties[name]["spatial"]:

                    # Selection of coordinate x, y, z for spatial symmetry
                    coordinate: int = spatial - 3 * int(spatial/3)
                    atom: int = int(spatial/3)
                    
                    if atom >= number_atoms:
                        raise ValueError(f"***Error \n\n\
                            atom {atom} doesn't exist") 
    
                    if coordinate == 0:
                        spatial_sym: int = 0
                    elif coordinate == 1:
                        spatial_sym: int = 1
                    elif coordinate == 2:
                        spatial_sym: int = 2
                    else:
                        raise ValueError("*** Error\n\n \
                            spatial sym doesn't exist")

                    for m_component in properties[name]["magnetic"]:

                        if type(spatial) == int:
                            integral_name: str = (name.lower() + " " +
                            str(spatial + 1) + " " + magnetic_components[m_component])
                            magnetic_xyz: int = m_component
                        else:
                            integral_name: str = (name.lower() + " " +
                            str(spatial + 1) + " "  + m_component)
                            magnetic_xyz: int = (list(magnetic_components.keys())
                            [list(magnetic_components.values()).index(m_component)])

                        symmetry[
                            integral_name
                        ] = integral_symmetry[name.lower()]

                        integrals[integral_name] = h1i(
                            charge = self._charge,
                            coord = self._coord,
                            exp = self._exp,
                            center = self._center,
                            lx = self._lx,
                            ly = self._ly,
                            lz = self._lz,
                            name = name,
                            output = output,
                            magnetic_xyz = magnetic_xyz,
                            spatial_sym = spatial_sym,
                            atom = atom,
                            gauge = properties[name]["gauge"]
                        )

        # Print integral
        if output > 10:
            for atomic_integrals_name in integrals.keys():
                print_triangle_matrix(
                    vector_to_matrix(
                        len(self._exp),
                        integrals[
                            atomic_integrals_name
                        ],
                        symmetry[
                            atomic_integrals_name
                        ],
                    ),
                    atomic_integrals_name,
                    symmetry[
                        atomic_integrals_name
                    ]
                )
        return integrals


if __name__ == "__main__":
    from wave_function import *

    wfn = wave_function("io/H2.molden")

    s = eint(wfn.build_wfn_array())

    s.integration(["ozke"],
                  #["overlap", "pot", "angmom"], 
                  {
                  "pot":{"atoms":[0, 1]}, 
                  "angmom":{"magnetic":[0, 1, 2], "gauge":[0.0, 0.0, 1.404552358700]},
                  "sd":{"spatial":[0,1,2,3,4,5], "magnetic":[0,1,2]},
                  "fc":{"atoms":[0,1]},
                  "nelfld":{"spatial":[0,1,2,3,4,5]},
                  "diplen":{"rdipole":[0.0,0.0,0.0],"magnetic":[0,1,2]},
                  "dipvel":{"magnetic":[0,1,2]},
                  "pso":{"spatial":[0,1,2,3,4,5]},
                  "nstcgo":{"spatial":[0,1,2,3,4,5],"magnetic":[0,1,2], "gauge":[0.0, 0.0, 1.404552358700]},
                  "dnske":{"spatial":[0,1,2,3,4,5],"magnetic":[0,1,2], "gauge":[0.0, 0.0, 1.404552358700]},
                  "psoke":{"spatial":[0,1,2,3,4,5]},
                  "psooz":{"spatial":[0,1,2,3,4,5],"magnetic":[0,1,2], "gauge":[0.0, 0.0, 1.404552358700]},
                  "ozke":{"magnetic":[0,1,2], "gauge":[0.0, 0.0, 1.404552358700]},
                  },
                  12)
