from typing import Union

from libr import *
from gradient_property import *
from coulomb_exchange import *
from principal_propagator import *
from lineal_response import *
from quadratic_response import *

multiplicity_string = {1: "Singlet", 3: "Triplet"}


class response:
    def __init__(
        self,
        wf: wave_function,
    ):
        """
        Reponse object

        This object drives response calculation like rpa, hrpa, etc.
        Build:
            *) Main Propagator
            *) Gradient Properties
            *) Lineal, Quadractic, ... response calculate, <<A; B>>

        Response calculation can be initialized from:
            * Wave Function Object, with the method reponse_wf
            * Using information neccesary for calulation: Molecular
            * Also, is possible combined methodologies
            Orbital Energies, Gradient Properties and Two-Body
            Integrals

        Gradient Properties values is a dict(list)
            Example:
                [
                    "Property 1": [a, b, c, ...], <-- Valors for property one
                    "Property 2": [d, e, f, ...]  <-- Valors for property two
                ]
        Args:
        ----
        wf (wave_function): Wave Function object
        moe (list): Molecular Orbital Energies sort in a list
        two_integrals (list or array): Molecular orbitals two-body integrals
        level (str): Approximation level pzoa, rpa, hrpa, ...
        properties (list(list(str))): Name of the properties involved
                in the response in lists
            Example:
                [
                    ["FC 1", "FC 2"], <-- Lineal Response
                    ...
                ]
        pp_multiplicity (list(list)): Multiplicity to build the principal propagator
        property_multiplicity (list(list)): Multiplicity for each property
        verbose (int): Print level
        """

        if not wf:  # or (not moe and not gradient_properties and not two_integrals
            # and not moe):
            raise ValueError(
                "** ERROR \n\n\
                    The information is not enogh to calculate. It is necessary\n\
                    wave function object or information by separated. In both\n\
                    cases is neccesary indicate the multiplicity and properties\n\
                    for each repsonse.\
                    "
            )

        self._wf: wave_function = wf

        self.principal_propagator = {"singlet": np.array([]), "triplet": np.array([])}
        self.q_principal_propagator = {"singlet": np.array([])}

    ################################################################################################
    # ATTRIBUTES
    ################################################################################################
    @property
    def principal_propagator_multiplicites(self) -> list:
        "Multiplicites"
        return [
            m if isinstance(m, int) else m.lower()
            for ms in self._pp_multiplicity
            for m in ms
        ]

    @property
    def principal_propagator_singlet(self) -> bool:
        "Singlet Multiplicity"
        if 1 in self.principal_propagator_multiplicites:
            return True
        if "singlet" in self.principal_propagator_multiplicites:
            return True
        return False

    @property
    def principal_propagator_triplet(self) -> bool:
        "Triplet Multiplicity"
        if 3 in self.principal_propagator_multiplicites:
            return True
        if "triplet" in self.principal_propagator_multiplicites:
            return True
        return False

    @property
    def property_multiplicity(self) -> dict:
        return {
            p: m if isinstance(m, str) else multiplicity_string[m]
            for pro, ms in zip(self._properties, self._p_multiplicity)
            for p, m in zip(pro, ms)
        }

    @property
    def lineal_property_multiplicity(self) -> list:
        return [
            m if isinstance(m, str) else multiplicity_string[m]
            for ms in self._p_multiplicity
            if len(ms) == 2
            for m in ms
        ]

    @property
    def quadratic_property_multiplicity(self) -> list:
        return [
            m if isinstance(m, str) else multiplicity_string[m]
            for ms in self._p_multiplicity
            if len(ms) == 3
            for m in ms
        ]

    @property
    def properties(self) -> list:
        "Properties different to calcualte"
        return [name for name in {name for p in self._properties for name in p}]

    @property
    def lineal_properties(self) -> list:
        "Properties different to calcualte to lineal response"
        return [name for p in self._properties if len(p) == 2 for name in p]

    @property
    def quadratic_properties(self) -> list:
        "Properties different to calcualte to lineal quadratic"
        return [name for p in self._properties if len(p) == 3 for name in p]

    @property
    def type_response(self) -> list:
        "Response Type"
        response = []
        for property in self._properties:
            if len(property) > 3:
                raise ValueError(
                    "***Error \n\n\
                        Response isn't implemete. The lineal and quadractic response\n\
                        is unique implemeted, i.e., reponse between 2 properties. Then ,\n\
                        the reponse array only might be [[A,B],[C,D,E],...]"
                )
            response.append(response_type[str(len(property))])
        return response

    @property
    def lineal_multiplicity_string(self) -> list:
        ms = []
        for multiplicity in self._pp_multiplicity:
            if len(multiplicity) == 1:
                if isinstance(multiplicity[0], int):
                    ms.append(multiplicity_string[multiplicity[0]])
                else:
                    ms.append(multiplicity[0])
        return ms

    @property
    def quadratic_multiplicity_string(self) -> list:
        ms = []
        for multiplicity in self._pp_multiplicity:
            if len(multiplicity) == 3:
                for m in multiplicity:
                    if isinstance(m, int):
                        ms.append(multiplicity_string[m])
                    else:
                        ms.append(m)
        return ms

    ################################################################################################
    # METHODS
    ################################################################################################
    def rpa(self, io: scratch, driver_time: drv_time):
        """
        Calculate reponse at random phase approximation

        Args:
        ----
            io (object:scratch): Driver to driver the output and binary files
            time_object (drv_time): Manage time calculation
        """

        io.write_output(information="RPA", type=1, title_type=1)
        start: float = time()

        countl = countlp = 0
        countq = countqp = 0
        for iresponse, property in enumerate(self._properties):

            # Type of response
            io.write_output("Response Calculation Specification")
            io.write_output(
                f"     Response Type: {self.type_response[iresponse].title()}"
            )
            # Properties pair
            operator_a: list = []
            operator_b: list = []
            operator_c: list = []
            for name in self.gpvs.keys():
                if property[0] in name:
                    operator_a.append(name)
                if property[1] in name:
                    operator_b.append(name)
                if self.type_response[iresponse] == "quadratic" and property[2] in name:
                    operator_c.append(name)

            if self.type_response[iresponse] == "lineal":
                properties: str = (
                    self.lineal_properties[countlp]
                    + ", "
                    + self.lineal_properties[countlp + 1]
                )
                io.write_output(f"     Properties: {properties}")
                multiplicites: str = (
                    self.lineal_property_multiplicity[countlp]
                    + ", "
                    + self.lineal_property_multiplicity[countlp + 1]
                )
                io.write_output(f"     Property Multiplicity: {multiplicites}")
                countlp += 2
                io.write_output(
                    f"     Principal Propagator Multiplicity: {self.lineal_multiplicity_string[countl]}"
                )
                for index_a, name_a in enumerate(operator_a):
                    for index_b, name_b in enumerate(operator_b):
                        if index_a > index_b and name_a.split()[0] == name_b.split()[0]:
                            continue
                        io.write_output(
                            f"     Lineal Response: <<{name_a.upper()};{name_b.upper()}>>"
                        )

                # Lineal Response calculate
                responses_values = calculate_lineal_reponse(
                    io=io,
                    time_object=driver_time,
                    operator_a=operator_a,
                    operator_b=operator_b,
                    n_mo_occ=self._wf.mo_occ,
                    n_mo_virt=self._wf.mo_virt,
                    pp_multiplicity=self.lineal_multiplicity_string[countl].lower(),
                    gpvs=self.gpvs,
                    verbose=self._verbose,
                )
                countl += 1

            elif self.type_response[iresponse] == "quadratic":
                #
                properties = (
                    self.quadratic_properties[countqp]
                    + ", "
                    + self.quadratic_properties[countqp + 1]
                    + ", "
                    + self.quadratic_properties[countqp + 2]
                )
                io.write_output(f"     Properties: {properties}")
                multiplicites = (
                    self.quadratic_property_multiplicity[countqp]
                    + ", "
                    + self.quadratic_property_multiplicity[countqp + 1]
                    + ", "
                    + self.quadratic_property_multiplicity[countqp + 2]
                )
                io.write_output(f"     Property Multiplicity: {multiplicites}")
                countqp += 3
                multiplicites = (
                    self.quadratic_multiplicity_string[countq]
                    + ", "
                    + self.quadratic_multiplicity_string[countq + 1]
                    + ", "
                    + self.quadratic_multiplicity_string[countq + 2]
                )
                io.write_output(
                    f"     Principal Propagator Multiplicity: {multiplicites}"
                )
                for index_a, name_a in enumerate(operator_a):
                    for index_b, name_b in enumerate(operator_b):
                        if index_a > index_b and name_a.split()[0] == name_b.split()[0]:
                            continue
                        for index_c, name_c in enumerate(operator_c):
                            if (
                                index_b > index_c
                                and name_b.split()[0] == name_b.split()[0]
                            ):
                                continue
                            io.write_output(
                                f"     Quadratic Response: <<{name_a.upper()};{name_b.upper()},{name_c.upper()}>> \n"
                            )
                #
                responses_values = calculate_quadratic_response(
                    io=io,
                    time_object=driver_time,
                    operator_a=operator_a,
                    operator_b=operator_b,
                    operator_c=operator_c,
                    n_mo_occ=self._wf.mo_occ,
                    n_mo_virt=self._wf.mo_virt,
                    pp_multiplicity_a=self.quadratic_multiplicity_string[
                        countq
                    ].lower(),
                    pp_multiplicity_b=self.quadratic_multiplicity_string[
                        countq + 1
                    ].lower(),
                    pp_multiplicity_c=self.quadratic_multiplicity_string[
                        countq + 2
                    ].lower(),
                    mo_occupied=self.mo_occupied,
                    mo_virtuals=self.mo_virtuals,
                    gpvs=self.gpvs,
                    verbose=self._verbose,
                )
                countq += 3

        driver_time.add_name_delta_time(name="RPA", delta_time=(time() - start))
        return responses_values

    #################################
    def drv_reponse_calculation(
        self,
        principal_propagator_approximation: str = "rpa",
        properties: list = [],
        property_multiplicity: list = [],
        pp_multiplicity: list = [],
        gauge: list = [0.0, 0.0, 0.0],
        verbose: int = 0,
        verbose_integrals: int = -1,
    ) -> dict[str, float]:
        """
        Manage of response calculation
        Args:
            principal_propagator_approximation (str): Approximation to principal propagator: PZOA, RPA, HRPA
                                                    by default to do RPA calculation
            gauge  (list): Gauge coordinate
            verbise (int): Print level
                    < 10: minimun
                    > 10: Time with more details
                    > 30: Gradient Vector properties, ocuppied molecular integrals obirtals, virtuals molecular integrals orbitals
                    > 50: Exchange and coulomb integrals, principal propagator and its inverse
            verbose_integrals (int): Print level to hermite module
        """

        ## Instance external objects
        # - Scratch
        io = self._wf._driver_scratch
        # - Driver Time
        driver_time = self._wf._driver_time
        ##
        start: float = time()

        io.write_output(information=f"REPONSE MODULE", type=1)

        self._properties = properties
        self._p_multiplicity: list = property_multiplicity
        if len(self._p_multiplicity) == 0:
            self._p_multiplicity = []
            for property in self._properties:
                temp: list = []
                for p in property:
                    temp.append(property_multiplicities[p.lower().split()[0]])
                self._p_multiplicity.append(temp)

        self._pp_multiplicity: list = pp_multiplicity
        if len(self._pp_multiplicity) == 0:
            self._pp_multiplicity = []
            for pm in self._p_multiplicity:
                if len(pm) == 2:
                    self._pp_multiplicity.append([pm[1]])
                elif len(pm) == 3:
                    self._pp_multiplicity.append(pm[0:3])
                else:
                    raise ValueError(
                        "***ERROR\n\n Response is not implemeted.\n\n\
                        Only is implemeted lineal and quadratic response."
                    )

        self._verbose = verbose

        # - Molecular Orbital Energies
        self._moe = np.array(self._wf.mo_energies)

        # - Two--Body Atomic Orbitalas
        calculate_integral = eint(self._wf)
        if not io._hermite_ao2b_binary.exists():
            calculate_integral.integration_twobody(
                integrals_names=["e2pot"],
                verbose=verbose_integrals,
            )
        # atomic one-body integrals
        calculate_integral.integration_onebody(
            integrals_names=self.properties, verbose=verbose_integrals, gauge=gauge
        )

        if "quadratic" in self.type_response:
            quadratic: bool = True
        else:
            quadratic = False

        # Gradient Property Vector and Average Values
        list_1b_integrals: list = calculate_integral.list_1b_integrals
        self.mo_occupied, self.mo_virtuals, self.gpvs = drv_gradient_property_vector(
            io=io,
            driver_time=driver_time,
            wf=self._wf,
            list_1b_integrals=list_1b_integrals,
            quadratic=quadratic,
            verbose=self._verbose,
        )

        # Principal Propagator
        for ms in self.lineal_multiplicity_string:
            if not ms.lower() in [1, 3, "singlet", "triplet"]:
                raise ValueError(
                    f"***ERROR\n\n\
                    The multiplicity {ms} is not implemented for principal propagator"
                )

        # Lineal principal propagator
        singlet_activate: bool = False
        triplet_activate: bool = False
        if not io._principal_propagator.exists():
            for ms in self._pp_multiplicity:
                if 1 in ms:
                    singlet_activate = True
                if 3 in ms:
                    triplet_activate = True
                if singlet_activate and triplet_activate:
                    break
        elif not io.binary(
            file=io._principal_propagator, io="f", label="singlet"
        ) or not io.binary(file=io._principal_propagator, io="f", label="triplet"):
            for ms in self._pp_multiplicity:
                if 1 in ms and not io.binary(
                    file=io._principal_propagator, io="f", label="singlet"
                ):
                    singlet_activate = True
                if 3 in ms and not io.binary(
                    file=io._principal_propagator, io="f", label="triplet"
                ):
                    triplet_activate = True
                if singlet_activate and triplet_activate:
                    break

        if singlet_activate or triplet_activate:
            # - Coulomb and Exchange
            if not io._exchange_coulomb.exists():
                get_coulomb_exchange_integrals(
                    io=io,
                    driver_time=driver_time,
                    wf=self._wf,
                    verbose=self._verbose,
                )
            drv_principal_propagator(
                io=io,
                driver_time=driver_time,
                n_mo_occ=self._wf.mo_occ,
                n_mo_virt=self._wf.mo_virt,
                moe=self._moe,
                multiplicity_pp={
                    "singlet": singlet_activate,
                    "triplet": triplet_activate,
                },
                tp_inv=0,
                verbose=self._verbose,
            )

        # Run Response
        if principal_propagator_approximation.lower() == "rpa":
            responses_values = self.rpa(
                io=io,
                driver_time=driver_time,
            )

        driver_time.add_name_delta_time(
            name="Response Calculation", delta_time=(time() - start)
        )
        io.write_output(type=2, drv_time=driver_time)
        driver_time.reset

        io.write_output(information=f"END REPONSE CALCULATION", type=1)

        return responses_values


if __name__ == "__main__":
    wfn = wave_function(
        "../tests/molden_file/LiH.molden",
        scratch_path="/home1/scratch",
        #job_folder="160922134451",
    )
    r = response(wfn)
    # r.drv_reponse_calculation(principal_propagator_approximation="rpa",
    #         properties = [["angmom x","fc 1","spinorbit x"],["angmom x","sd 1 x","spinorbit x"],["angmom x","sd 1 z","spinorbit z"],["angmom y","sd 2 y","spinorbit y"]
    #         ,["angmom z","sd 3 x","spinorbit x"],["angmom x", "pso 1", "massvelo"],["angmom x", "pso 1", "darwin"]],
    #                             pp_multiplicity=[[1,3,3],[1,3,3],[1,3,3],[1,3,3],[1,3,3],[1,1,1],[1,1,1]],gauge=[0.000,0.0000,0.0586476414],
    #                             verbose=11)

    run = True
    if run:
        r.drv_reponse_calculation(
            principal_propagator_approximation="rpa",
            properties=[
                ["fc", "kinetic"],
                ["fc 1", "laplacian xx"],
                ["fc 1", "laplacian yy"],
                ["fc 1", "laplacian zz"],
                ["angmom x", "fc 1", "spinorbit x"],
                ["angmom x", "fc 2", "spinorbit x"],
                ["fc", "fc"],
                ["spinorbit", "fc"],
                ["spinorbit", "spinorbit"],
                ["pso", "pso"],
            ],
            # gauge=[0.0,0.0,1.4045523587],
            # gauge = [0.000, 0.0000, -0.545857052], #Li pople
            verbose=11,
            verbose_integrals=1,
        )
    else:
        a = 0
        r.drv_reponse_calculation(
            principal_propagator_approximation="rpa",
            properties=[
                ["angmom x", "pso " + str(1 + 3 * a)],
                ["angmom y", "pso " + str(2 + 3 * a)],
                ["angmom z", "pso " + str(3 + 3 * a)],
                ["pso 1", "ozke x"],
                ["pso 2", "ozke y"],
                ["pso 3", "ozke z"],
                ["angmom x", "psoke 1"],
                ["angmom y", "psoke 2"],
                ["angmom z", "psoke 3"],
                ["laplacian xx", "sd " + str(1 + 3 * a) + " x"],
                ["laplacian xy", "sd " + str(1 + 3 * a) + " x"],
                ["laplacian xz", "sd " + str(1 + 3 * a) + " x"],
                ["laplacian xx", "sd " + str(1 + 3 * a) + " y"],
                ["laplacian xy", "sd " + str(1 + 3 * a) + " y"],
                ["laplacian xz", "sd " + str(1 + 3 * a) + " y"],
                ["laplacian xx", "sd " + str(1 + 3 * a) + " z"],
                ["laplacian xy", "sd " + str(1 + 3 * a) + " z"],
                ["laplacian xz", "sd " + str(1 + 3 * a) + " z"],
                ["laplacian xy", "sd " + str(2 + 3 * a) + " x"],
                ["laplacian yy", "sd " + str(2 + 3 * a) + " x"],
                ["laplacian yz", "sd " + str(2 + 3 * a) + " x"],
                ["laplacian xy", "sd " + str(2 + 3 * a) + " y"],
                ["laplacian yy", "sd " + str(2 + 3 * a) + " y"],
                ["laplacian yz", "sd " + str(2 + 3 * a) + " y"],
                ["laplacian xy", "sd " + str(2 + 3 * a) + " z"],
                ["laplacian yy", "sd " + str(2 + 3 * a) + " z"],
                ["laplacian yz", "sd " + str(2 + 3 * a) + " z"],
                ["laplacian xz", "sd " + str(3 + 3 * a) + " x"],
                ["laplacian yz", "sd " + str(3 + 3 * a) + " x"],
                ["laplacian zz", "sd " + str(3 + 3 * a) + " x"],
                ["laplacian xz", "sd " + str(3 + 3 * a) + " y"],
                ["laplacian yz", "sd " + str(3 + 3 * a) + " y"],
                ["laplacian zz", "sd " + str(3 + 3 * a) + " y"],
                ["laplacian xz", "sd " + str(3 + 3 * a) + " z"],
                ["laplacian yz", "sd " + str(3 + 3 * a) + " z"],
                ["laplacian zz", "sd " + str(3 + 3 * a) + " z"],
                ["kinetic", "fc 1"],
                ["fc " + str(1 + a), "sofiel xx"],
                ["fc " + str(1 + a), "sofiel yy"],
                ["fc " + str(1 + a), "sofiel zz"],
                ["sd " + str(1 + 3 * a) + " x", "sofiel xx"],
                ["sd " + str(1 + 3 * a) + " x", "sofiel yx"],
                ["sd " + str(1 + 3 * a) + " x", "sofiel zx"],
                ["sd " + str(1 + 3 * a) + " y", "sofiel xx"],
                ["sd " + str(1 + 3 * a) + " y", "sofiel yx"],
                ["sd " + str(1 + 3 * a) + " y", "sofiel zx"],
                ["sd " + str(1 + 3 * a) + " z", "sofiel xx"],
                ["sd " + str(1 + 3 * a) + " z", "sofiel yx"],
                ["sd " + str(1 + 3 * a) + " z", "sofiel zx"],
                ["sd " + str(2 + 3 * a) + " x", "sofiel xy"],
                ["sd " + str(2 + 3 * a) + " x", "sofiel yy"],
                ["sd " + str(2 + 3 * a) + " x", "sofiel zy"],
                ["sd " + str(2 + 3 * a) + " y", "sofiel xy"],
                ["sd " + str(2 + 3 * a) + " y", "sofiel yy"],
                ["sd " + str(2 + 3 * a) + " y", "sofiel zy"],
                ["sd " + str(2 + 3 * a) + " z", "sofiel xy"],
                ["sd " + str(2 + 3 * a) + " z", "sofiel yy"],
                ["sd " + str(2 + 3 * a) + " z", "sofiel zy"],
                ["sd " + str(3 + 3 * a) + " x", "sofiel xz"],
                ["sd " + str(3 + 3 * a) + " x", "sofiel yz"],
                ["sd " + str(3 + 3 * a) + " x", "sofiel zz"],
                ["sd " + str(3 + 3 * a) + " y", "sofiel xz"],
                ["sd " + str(3 + 3 * a) + " y", "sofiel yz"],
                ["sd " + str(3 + 3 * a) + " y", "sofiel zz"],
                ["sd " + str(3 + 3 * a) + " z", "sofiel xz"],
                ["sd " + str(3 + 3 * a) + " z", "sofiel yz"],
                ["sd " + str(3 + 3 * a) + " z", "sofiel zz"],
                ["angmom x", "fc " + str(1 + a), "spinorbit x"],
                ["angmom y", "fc " + str(1 + a), "spinorbit y"],
                ["angmom z", "fc " + str(1 + a), "spinorbit z"],
                ["angmom x", "sd " + str(1 + 3 * a) + " x", "spinorbit x"],
                ["angmom x", "sd " + str(1 + 3 * a) + " y", "spinorbit y"],
                ["angmom x", "sd " + str(1 + 3 * a) + " z", "spinorbit z"],
                ["angmom y", "sd " + str(2 + 3 * a) + " x", "spinorbit x"],
                ["angmom y", "sd " + str(2 + 3 * a) + " y", "spinorbit y"],
                ["angmom y", "sd " + str(2 + 3 * a) + " z", "spinorbit z"],
                ["angmom z", "sd " + str(3 + 3 * a) + " x", "spinorbit x"],
                ["angmom z", "sd " + str(3 + 3 * a) + " y", "spinorbit y"],
                ["angmom z", "sd " + str(3 + 3 * a) + " z", "spinorbit z"],
                ["angmom x", "pso " + str(1 + 3 * a), "massvelo"],
                ["angmom y", "pso " + str(2 + 3 * a), "massvelo"],
                ["angmom z", "pso " + str(3 + 3 * a), "massvelo"],
                ["angmom x", "pso " + str(1 + 3 * a), "darwin"],
                ["angmom y", "pso " + str(2 + 3 * a), "darwin"],
                ["angmom z", "pso " + str(3 + 3 * a), "darwin"],
                ["nstcgo " + str(1 + 3 * a) + " x", "darwin"],
                ["nstcgo " + str(2 + 3 * a) + " y", "darwin"],
                ["nstcgo " + str(3 + 3 * a) + " z", "darwin"],
                ["nstcgo " + str(1 + 3 * a) + " x", "massvelo"],
                ["nstcgo " + str(2 + 3 * a) + " y", "massvelo"],
                ["nstcgo " + str(3 + 3 * a) + " z", "massvelo"],
            ],
            gauge=[0.000, 0.0000, -0.545857052],
            verbose=11,
        )

    # r.drv_reponse_calculation(principal_propagator_approximation="rpa", properties = [["kinetic","kinetic","kinetic"]],
    #                             pp_multiplicity=[[3,3,3]],
    #                             verbose=11)
