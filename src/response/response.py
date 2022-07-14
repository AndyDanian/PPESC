

from libr import *
from gradient_property import *
from coulomb_exchange import *
from principal_propagator import *
from lineal_response import *
from quadratic_response import *

multiplicity_string = {1:"Singlet", 3:"Triplet"}

class response():
    def __init__(self, wf: wave_function = None, moe: list = None,
                    gradient_properties: dict = None,
                    coulomb_integrals: np.array = None,
                    exchange_integrals: np.array = None):
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

        Args:
        ----
        wf (wave_function): Wave Function object
        moe (list): Molecular Orbital Energies sort in a list
        gradient_properties (dict(list)): Gradient Properties values in a list
            Example:
                [
                    "Property 1": [a, b, c, ...], <-- Valors for property one
                    "Property 2": [d, e, f, ...]  <-- Valors for property two
                ]
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

        if not wf: # or (not moe and not gradient_properties and not two_integrals
            #and not moe):
            raise ValueError("** ERROR \n\n\
                    The information is not enogh to calculate. It is necessary\n\
                    wave function object or information by separated. In both\n\
                    cases is neccesary indicate the multiplicity and properties\n\
                    for each repsonse.\
                    ")

        self._wf = wf
        self._moe = moe
        self._gp = gradient_properties

        if not coulomb_integrals: self._coulomb_integrals = np.array([])
        else:   self._coulomb_integrals = coulomb_integrals
        if not exchange_integrals: self._exchange_integrals = np.array([])
        else: self._exchange_integrals = exchange_integrals

        self.principal_propagator = {"singlet": np.array([]), "triplet": np.array([])}

        if not self._gp:
            self._gp = {}
    ################################################################################################
    # ATTRIBUTES
    ################################################################################################
    @property
    def principal_propagator_multiplicites(self) -> list:
        "Multiplicites"
        return [m if isinstance(m, int) else m.lower() for ms in self._pp_multiplicity for m in ms]
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
        return {p: m if isinstance(m, str) else multiplicity_string[m]
                for pro, ms in zip(self._properties, self._p_multiplicity) for p, m in zip(pro, ms)}
    @property
    def lineal_property_multiplicity(self) -> dict:
        return [m if isinstance(m, str) else multiplicity_string[m]
                for ms in self._p_multiplicity
                if len(ms) == 2 for m in ms]
    @property
    def quadratic_property_multiplicity(self) -> dict:
        return [m if isinstance(m, str) else multiplicity_string[m]
                for ms in self._p_multiplicity
                if len(ms) == 3 for m in ms]
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
                raise ValueError("***Error \n\n\
                        Response isn't implemete. The lineal and quadractic response\n\
                        is unique implemeted, i.e., reponse between 2 properties. Then ,\n\
                        the reponse array only might be [[A,B],[C,D,E],...]")
            response.append(response_type[str(len(property))])
        return response
    @property
    def lineal_multiplicity_string(self)-> list:
        ms = []
        for multiplicity in self._pp_multiplicity:
            if len(multiplicity) == 1:
                if isinstance(multiplicity[0], int):
                    ms.append(multiplicity_string[multiplicity[0]])
                else:
                    ms.append(multiplicity[0])
        return ms
    @property
    def quadratic_multiplicity_string(self)-> list:
        ms = []
        for multiplicity in self._pp_multiplicity:
            if len(multiplicity) == 2:
                for m in multiplicity:
                    if isinstance(m, int):
                        ms.append(multiplicity_string[m])
                    else:
                        ms.append(m)
        return ms
    ################################################################################################
    # METHODS
    ################################################################################################
    def rpa(self, driver_time: drv_time = None):
        """
        Calculate reponse at random phase approximation

        Args:
        ----
        verbose_integrals (int): Print level for integrals calculation
        """

        print_subtitle(name = "RPA")

        if self._verbose > 10:
            start = time()

        countl = countlp = 0
        countq = countqp = 0
        for iresponse, property in enumerate(self._properties):

            # Type of response
            print_ljust(name = "Response Calculation Specification")
            print(f"     Response Type: {self.type_response[iresponse].title()}")
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
                print(f"     Properties: {self.lineal_properties[countlp]}, {self.lineal_properties[countlp + 1]}")
                print("     Property Multiplicity: ",self.lineal_property_multiplicity[countlp],", ",
                                                    self.lineal_property_multiplicity[countlp + 1])
                countlp += 2
                print(f"     Principal Propagator Multiplicity: {self.lineal_multiplicity_string[countl]}")
                for index_a, name_a in enumerate(operator_a):
                    for index_b, name_b in enumerate(operator_b):
                        if index_a > index_b and name_a.split()[0] == name_b.split()[0]:
                            continue
                        print(f"     Lineal Response: <<{name_a.upper()};{name_b.upper()}>>")
                #Response calculate
                responses_values = calculate_lineal_reponse(
                operator_a = operator_a, operator_b = operator_b,
                n_mo_occ = self._wf.mo_occ, n_mo_virt = self._wf.mo_virt,
                principal_propagator = self.principal_propagator[self.lineal_multiplicity_string[countl].lower()],
                gpvs = self.gpvs, time_object = driver_time, verbose = self._verbose)
                countl += 1

            elif self.type_response[iresponse] == "quadratic":
            #
                print(f"     Properties: {self.quadratic_properties[countqp]}, {self.quadratic_properties[countqp + 1]} , {self.quadratic_properties[countqp + 2]}")
                print("     Property Multiplicity: ",self.quadratic_property_multiplicity[countqp],", ",
                                                    self.quadratic_property_multiplicity[countqp + 1],", ",
                                                    self.quadratic_property_multiplicity[countqp + 2])
                countqp += 3
                print("     Principal Propagator Multiplicity: ",self.quadratic_multiplicity_string[countq],", ",
                                                                self.quadratic_multiplicity_string[countq + 1])
                for index_a, name_a in enumerate(operator_a):
                    for index_b, name_b in enumerate(operator_b):
                        if index_a > index_b and name_a.split()[0] == name_b.split()[0]:
                            continue
                        for index_c, name_c in enumerate(operator_c):
                            if index_b > index_c and name_b.split()[0] == name_b.split()[0]:
                                continue
                            print(f"     Quadratic Response: <<{name_a.upper()};{name_b.upper()},{name_c.upper()}>>")
                #
                responses_values = calculate_quadratic_response(
                operator_a = operator_a, operator_b = operator_b, operator_c = operator_c,
                n_mo_occ = self._wf.mo_occ, n_mo_virt = self._wf.mo_virt,
                principal_propagator_a = self.principal_propagator[self.quadratic_multiplicity_string[countq].lower()],
                principal_propagator_b = self.principal_propagator[self.quadratic_multiplicity_string[countq + 1].lower()],
                avs = self.avs, mo_occupied = self.mo_occupied, mo_virtuals = self.mo_virtuals, gpvs = self.gpvs,
                time_object = driver_time, verbose = self._verbose)
                countq += 2

        if self._verbose > 10:
            driver_time.add_name_delta_time(name = "RPA", delta_time = (time() - start))

        return responses_values

    #################################
    def drv_reponse_calculation(self, principal_propagator_approximation: str = "rpa", properties: list = None,
                                property_multiplicity: list = None, pp_multiplicity: list = None,
                                gaugeo: list = None,
                                verbose: int = 0, verbose_integrals: int = -1):
        """
        Manage of response calculation
        Args:
            principal_propagator_approximation (str): Approximation to principal propagator: PZOA, RPA, HRPA
                                                    by default to do RPA calculation
            gaugeo (list): Gauge coordinate
            verbise (int): Print level
            verbose_integrals (int): Print level to hermite module
        """

        self._properties = properties
        self._p_multiplicity = property_multiplicity
        if self._p_multiplicity is None:
            self._p_multiplicity: list = []
            for property in self._properties:
                temp: list = []
                for p in property:
                    temp.append(property_multiplicities[p.lower().split()[0]])
                self._p_multiplicity.append(temp)

        self._pp_multiplicity = pp_multiplicity
        if self._pp_multiplicity is None:
            self._pp_multiplicity: list = []
            for pm in self._p_multiplicity:
                if len(pm) == 2:
                    self._pp_multiplicity.append([pm[1]])
                elif len(pm) == 3:
                    self._pp_multiplicity.append(pm[1:3])
                else:
                    raise ValueError("***ERROR\n\n Response is not implemeted.\n\n\
                        Only is implemeted lineal and quadratic response.")

        self._verbose = verbose

        if verbose >= 0:
            print_title(name = f"REPONSE CALCULATION")

        if self._verbose > 10:
            driver_time = drv_time()
            start = time()
        else:
            driver_time = None

        # - Molecular orbital energies
        if not self._moe:
            self._moe = self._wf.mo_energies

        # - Two integrals
        if not self._coulomb_integrals.size or not self._exchange_integrals.size:
            self._coulomb_integrals, self._exchange_integrals = get_coulomb_exchange_integrals(self._wf,
                                                time_object = driver_time,
                                                verbose = self._verbose,
                                                verbose_int = verbose_integrals)
        #Principal Propagator
        for ms in self.lineal_multiplicity_string:
            if not ms.lower() in [1,3,"singlet","triplet"]:
                raise ValueError(f"***ERROR\n\n\
                    The multiplicity {ms} is not implemented for principal propagator")

        if (not self.principal_propagator["singlet"].size
            or not self.principal_propagator["triplet"].size):
            dict_principal_propagator = drv_principal_propagator(driver_time = driver_time, moe = self._moe,
                                                        n_mo_occ = self._wf.mo_occ, n_mo_virt = self._wf.mo_virt,
                                                        coulomb = self._coulomb_integrals, exchange = self._exchange_integrals,
                                                        multiplicity_pp = {"singlet": self.principal_propagator_singlet,
                                                        "triplet": self.principal_propagator_triplet},
                                                        tp_inv = 0, verbose = self._verbose)
            self.principal_propagator.update(dict_principal_propagator)
        # Gradient Property Vector and Average Values
        if "quadratic" in self.type_response:
            average = True
        else:
            average = False
        self.avs, self.mo_occupied, self.mo_virtuals, self.gpvs = drv_gradient_property_vector(wf = self._wf,
                                            properties =self.properties,
                                            gpv_in = self._gp, driver_time = driver_time, gaugeo = gaugeo,
                                            properties_multiplicity = self.property_multiplicity, average = average,
                                            verbose = self._verbose, verbose_integrals = verbose_integrals)
        # Run Response
        if principal_propagator_approximation.lower() == "rpa":
            responses_values = self.rpa(driver_time = driver_time)

        if self._verbose > 10:
            driver_time.add_name_delta_time(name = "Response Calculation", delta_time = (time() - start))
            driver_time.printing()

        if verbose >= 0:
            print_title(name = f"END REPONSE CALCULATION")

        return responses_values

if __name__ == "__main__":
    wfn = wave_function("../tests/molden_file/LiH_s.molden")
    r = response(wfn)
    r.drv_reponse_calculation(principal_propagator_approximation="rpa", properties = [["kinetic","fc 1","fc 2"]], property_multiplicity=[[3,3,3]],
                                pp_multiplicity=[[3,3]],
                                verbose=31)

    # r.drv_reponse_calculation(principal_propagator_approximation="rpa", properties = [["kinetic","fc 1","fc 2"]], property_multiplicity=[[3,3,3]],
    #                             verbose=31)