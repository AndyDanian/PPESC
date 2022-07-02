

from gradient_property import *
from coulomb_exchange import *
from principal_propagator import *
from lineal_response import *
from quadratic_response import *
from libr import *

multiplicity_string = {1:"Singlet", 3:"Triplet"}

class response():
    def __init__(self, wf: wave_function = None, moe: list = None,
                    gradient_properties: dict = None, two_integrals: list or array = None,
                    properties: list = None, property_multiplicity: list = None,
                    pp_multiplicity: list = None,
                    verbose: int = 0):
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

        if not wf or (not moe and not gradient_properties and not properties
            and not multiplicity):
            raise ValueError("** ERROR \n\n\
                    The information is not enogh to calculate. It is necessary\n\
                    wave function object or information by separated. In both\n\
                    cases is neccesary indicate the multiplicity and properties\n\
                    for each repsonse.\
                    ")

        self._wf = wf
        self._properties = properties
        self._moe = moe
        self._gp = gradient_properties
        self._pp_multiplicity = pp_multiplicity
        self._p_multiplicity = property_multiplicity
        self._tintegrals = two_integrals
        self._verbose = verbose

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
    def lineal_singlet(self) -> bool:
        "Singlet Multiplicity"
        if 1 in self.principal_propagator_multiplicites:
            return True
        if "singlet" in self.principal_propagator_multiplicites:
            return True
        return False
    @property
    def lineal_triplet(self) -> bool:
        "Triplet Multiplicity"
        if 3 in self.principal_propagator_multiplicites:
            return True
        if "triplet" in self.principal_propagator_multiplicites:
            return True
        return False
    @property
    def property_multiplicity(self) -> dict:
        return {p: m for pro, ms in zip(self._properties, self._p_multiplicity) for p, m in zip(pro, ms)}
    @property
    def properties(self) -> list:
        "Properties different to calcualte"
        return [name for name in {name for p in self._properties for name in p}]
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
    def rpa(self, verbose_integrals: int = -1, verbose_average: int = -1):
        """
        Calculate reponse at random phase approximation

        Args:
        ----
        verbose_integrals (int): Print level for integrals calculation
        """

        print_title(name = "RESPONSE CALCULATION")

        if self._verbose > 10:
            driver_time = drv_time()
            start = time()
        else:
            driver_time = None

        print("\nLevel approximation: RPA \n")

        tr = self.type_response

        # - Molecular orbital energies
        if not self._moe:
            self._moe = self._wf.mo_energies
        # - Two integrals
        if not self._tintegrals:
            coulomb_integrals, exchange_integrals = get_coulomb_exchange_integrals(self._wf,
                                                time_object = driver_time,
                                                verbose = self._verbose,
                                                verbose_int = verbose_integrals)
        else:
            raise ValueError("Falta implementar que obtenga las integrales desde 2i dadas")
        #Principal Propagator
        for ms in self.lineal_multiplicity_string:
            if not ms.lower() in [1,3,"singlet","triplet"]:
                raise ValueError(f"***ERROR\n\n\
                    The multiplicity {ms} is not implemented for principal propagator")
        principal_propagator = drv_principal_propagator(driver_time = driver_time, moe = self._moe,
                                                        n_mo_occ = self._wf.mo_occ, n_mo_virt = self._wf.mo_virt,
                                                        coulomb = coulomb_integrals, exchange = exchange_integrals,
                                                        multiplicity_pp = {"singlet": self.lineal_singlet, "triplet": self.lineal_triplet},
                                                        tp_inv = 0, verbose = self._verbose)
        # Gradient Property Vector and Average Values
        if "quadratic" in self.quadratic_multiplicity_string:
            average = True
        else:
            average = False
        avs, mo_virtuals, gpvs = drv_gradient_property_vector(wf = self._wf, properties =self.properties, gpv_in = self._gp, driver_time = driver_time,
                                            properties_multiplicity = self.property_multiplicity, average = average,
                                            verbose = self._verbose, verbose_integrals = verbose_integrals)
        countl = 0
        countq = 0
        for iresponse, property in enumerate(self._properties):

            # Type of response
            print(f"     Response Type: {tr[iresponse]}")
            # Properties pair
            operator_a: list = []
            operator_b: list = []
            operator_c: list = []
            for name in gpvs.keys():
                if property[0] in name:
                    operator_a.append(name)
                if property[1] in name:
                    operator_b.append(name)
                if tr[iresponse] == "quadratic" and property[2] in name:
                    operator_c.append(name)

            if tr[iresponse] == "lineal":
                print(f"     Multiplicity: {self.lineal_multiplicity_string[countl]}")
                for index_a, name_a in enumerate(operator_a):
                    for index_b, name_b in enumerate(operator_b):
                        if index_a > index_b and name_a.split()[0] == name_b.split()[0]:
                            continue
                        print(f"     Lineal Response: <<{name_a.upper()};{name_b.upper()}>>")
                #Response calculate
                calculate_lineal_reponse(
                operator_a = operator_a, operator_b = operator_b,
                n_mo_occ = self._wf.mo_occ, n_mo_virt = self._wf.mo_virt,
                principal_propagator = principal_propagator[self.lineal_multiplicity_string[countl].lower()], gpvs = gpvs,
                time_object = driver_time, verbose = self._verbose)
                countl += 1
            elif tr[iresponse] == "quadratic":
                print(f"     Multiplicity: {self.quadratic_multiplicity_string[countq]} {self.quadratic_multiplicity_string[countq + 1]}")
                for index_a, name_a in enumerate(operator_a):
                    for index_b, name_b in enumerate(operator_b):
                        if index_a > index_b and name_a.split()[0] == name_b.split()[0]:
                            continue
                        for index_c, name_c in enumerate(operator_c):
                            if index_b > index_c and name_a.split()[0] == name_b.split()[0]:
                                continue
                            print(f"     Quadratic Response: <<{name_a.upper()};{name_b.upper()},{name_c.upper()}>>")
                #
                calculate_quadratic_response(
                operator_a = operator_a, operator_b = operator_b, operator_c = operator_c,
                n_mo_occ = self._wf.mo_occ, n_mo_virt = self._wf.mo_virt,
                principal_propagator_a = principal_propagator[self.quadratic_multiplicity_string[countq].lower()],
                principal_propagator_b = principal_propagator[self.quadratic_multiplicity_string[countq + 1].lower()],
                avs = avs, mo_virtuals = mo_virtuals, gpvs = gpvs,
                time_object = driver_time, verbose = self._verbose)
                countq += 2

        if self._verbose > 10:
            driver_time.add_name_delta_time(name = "Response Calculation", delta_time = (time() - start))
            driver_time.printing()
        print_title(name = f"END REPONSE CALCULATION")

if __name__ == "__main__":
    wfn = wave_function("../tests/molden_file/H2_s.molden")
    r = response(wfn, properties = [["fc","fc","fc"], ["fc","kinetic"]], property_multiplicity = [[3,3,3],[3,1]], pp_multiplicity=[[3,3],[1]], verbose=0)
    r.rpa()