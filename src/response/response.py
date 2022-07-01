

from gradient_property import *
from coulomb_exchange import *
from principal_propagator import *
from lineal_response import *
from libr import *

multiplicity_string = {1:"Singlet", 3:"Triplet"}

class response():
    def __init__(self, wf: wave_function = None, moe: list = None,
                    gradient_properties: dict = None, two_integrals: list or array = None,
                    properties: list = None, multiplicity: list = None,
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
        multiplicity (list(list)): Multiplicity of each property
        verbose (int): Print level
        """

        print_title(name = "RESPONSE CALCULATION")

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
        self._multiplicity = multiplicity
        self._tintegrals = two_integrals
        self._verbose = verbose

        if not self._gp:
            self._gp = {}
    ################################################################################################
    # ATTRIBUTES
    ################################################################################################
    @property
    def lineal_multiplicites(self) -> list:
        "Multiplicites"
        return [ms if isinstance(ms, int) else ms.lower() for ms in self._multiplicity]
    @property
    def lineal_singlet(self) -> bool:
        "Singlet Multiplicity"
        if 1 in self.lineal_multiplicites:
            return True
        if "singlet" in self.lineal_multiplicites:
            return True
        return False
    @property
    def lineal_triplet(self) -> bool:
        "Triplet Multiplicity"
        if 3 in self.lineal_multiplicites:
            return True
        if "triplet" in self.lineal_multiplicites:
            return True
        return False
    @property
    def properties(self) -> list:
        "Properties different to calcualte"
        return [name for name in {name for p in self._properties for name in p}]
    @property
    def lineal_multiplicity_properties(self) -> dict:
        "Properties with its multiplicity"
        return {p: self._multiplicity[int(count/2)] for count, pro in enumerate(self._properties) for p in pro}
    ################################################################################################
    # METHODS
    ################################################################################################
    def rpa(self, verbose_integrals: int = -1):
        """
        Calculate reponse at random phase approximation

        Args:
        ----
        verbose_integrals (int): Print level for integrals calculation
        """
        if self._verbose > 10:
            driver_time = drv_time()
            start = time()
        else:
            driver_time = None

        print("\nLevel approximation: RPA \n")

        #Principal Propagator
        # - Two integrals
        if not self._tintegrals:
            coulomb_integrals, exchange_integrals = get_coulomb_exchange_integrals(self._wf,
                                                    time_object = driver_time,
                                                    verbose = self._verbose,
                                                    verbose_int = verbose_integrals)
        else:
            raise ValueError("Falta implementar que obtenga las integrales desde 2i dadas")

        # - Molecular orbital energies
        if not self._moe:
            self._moe = self._wf.mo_energies

        # - Build PP
        multiplicity = {"singlet": self.lineal_singlet, "triplet": self.lineal_triplet}

        for ms in self.lineal_multiplicites:
            if not ms in [1,3,"singlet","triplet"]:
                raise ValueError(f"***ERROR\n\n\
                    The multiplicity is not implemented {ms}")

        principal_propagator = {}
        for name, ms in multiplicity.items():
            if ms:
                principal_propagator[name] = get_principal_propagator_lineal_rpa(
                                                            n_mo_occ = self._wf.mo_occ, n_mo_virt = self._wf.mo_virt,
                                                            moe = np.array(self._moe), coulomb = coulomb_integrals,
                                                            exchange = exchange_integrals,
                                                            multiplicity = name, tp_inv = 0,
                                                            time_object = driver_time,
                                                            verbose = self._verbose)

        # Gradient Property Vector
        gpvs: dict = {}
        for property in self.properties:
            if not self._gp or property not in self._gp.keys():
                temp_gpvs = gradient_property_vector_rpa(wf = self._wf,
                                    property = property, time_object = driver_time,
                                    verbose = self._verbose, verbose_int = verbose_integrals,
                                    multiplicity=self.lineal_multiplicity_properties[property])

                for name, value in temp_gpvs.items():
                    gpvs[name] = value
            else:
                gpvs = read_gradient_property_vector_rpa(wf = self._gp, property = property,
                                    verbose = self._verbose,
                                    multiplicity=self.lineal_multiplicity_properties[property])

        #Calculate of Response
        for iresponse, property in enumerate(self._properties):

            # Type of response
            if response_type[str(len(property))] != "lineal":
                raise ValueError("***Error \n\n\
                        Response isn't implemete. The lineal response is unique\n\
                        implemeted, i.e., reponse between 2 properties. Then, the\n\
                        reponse array only might be [[A,B],[C,D],...]")

            print(f"     Response Type: {response_type[str(len(property))]}")

            # Properties pair
            operator_a: list = []
            operator_b: list = []
            for name in gpvs.keys():
                if property[0] in name:
                    operator_a.append(name)
                if property[1] in name:
                    operator_b.append(name)

            if response_type[str(len(property))] == "lineal":
                for index_a, name_a in enumerate(operator_a):
                    for index_b, name_b in enumerate(operator_b):
                        if index_a > index_b and name_a.split()[0] == name_b.split()[0]:
                            continue
                        print(f"     Lineal Response: <<{name_a.upper()};{name_b.upper()}>>")

            #Multiplicity
            if not self._multiplicity[iresponse]:
                raise ValueError(f"*** ERROR\n\n\
                    Multiplicity array, {len(self._multiplicity)}, is lower than\n\
                    properties array, {len(self._properties)}.\n\
                    ")
            else:
                multiplicity: str = self._multiplicity[iresponse]

            if isinstance(multiplicity, int):
                multiplicity_s = multiplicity_string[multiplicity]
            else:
                multiplicity_s = multiplicity
            print(f"     Multiplicity: {multiplicity_s}")

            #Response calculate
            calculate_lineal_reponse(
            operator_a = operator_a, operator_b = operator_b,
            n_mo_occ = self._wf.mo_occ, n_mo_virt = self._wf.mo_virt,
            principal_propagator = principal_propagator[multiplicity_s.lower()], gpvs = gpvs,
            time_object = driver_time, verbose = self._verbose)

        if self._verbose > 10:
            driver_time.add_name_delta_time(name = "Response Calculation", delta_time = (time() - start))
            driver_time.printing()
        print_title(name = f"END REPONSE CALCULATION")

if __name__ == "__main__":
    wfn = wave_function("../tests/molden_file/H2_s.molden")
    r = response(wfn, properties = [["fc","fc"], ["fc","kinetic"]], multiplicity=[3,3], verbose=0)
    r.rpa()