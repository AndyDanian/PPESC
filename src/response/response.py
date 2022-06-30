

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
                    all_responses: bool = None, verbose: int = 0):
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
        all_response (bool): Activate the different response among one set
                            of operators
                            Example, if is True:
                            [A,B] then is get <<A;A>>, <<A;B>>, and <<B,B>>

                            Example, if is False:
                            [A,B] then is get <<A;B>>
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
        self._allr = all_responses
        self._verbose = verbose

        if not self._gp:
            self._gp = {}

    ################################################################################################
    # METHODS
    ################################################################################################
    def rpa(self, verbose_integrals: int = 0):
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
        for iresponse, response in enumerate(self._properties):

            if response_type[str(len(response))] != "lineal":
                raise ValueError("***Error \n\n\
                        Response isn't implemete. The lineal response is unique\n\
                        implemeted, i.e., reponse between 2 properties. Then, the\n\
                        reponse array only might be [[A,B],[C,D],...]")

            print(f"     Response Type: {response_type[str(len(response))]}")
            if response_type[str(len(response))] == "lineal":
                print(f"     Properties: <<{response[0].upper()};{response[1].upper()}>>")

            #Multiplicity
            if isinstance(self._multiplicity, list):
                if not self._multiplicity[iresponse]:
                    raise ValueError(f"*** ERROR\n\n\
                        Multiplicity array, {len(self._multiplicity)}, is lower than\n\
                        properties array, {len(self._properties)}.\n\
                        ")
                else:
                    multiplicity: list = self._multiplicity[iresponse]

            multiplicity_s = []
            for mult in multiplicity:
                if isinstance(mult, int):
                    multiplicity_s.append(multiplicity_string[mult])
                else:
                    multiplicity_s.append(mult)
            print(f"     Multiplicity: {multiplicity_s}")

            # Gradient Property Vector
            for count, property in enumerate(response):
                if count > 0 and property in response[:count]:
                    continue
                if not self._gp or property not in self._gp.keys():
                    temp_all_responses, gpvs = gradient_property_vector_rpa(wf = self._wf,
                                        property = property, time_object = driver_time,
                                        verbose = self._verbose, verbose_int = verbose_integrals,
                                        multiplicity=multiplicity[count])
                    if self._allr == None:
                        self._allr = temp_all_responses
                else:
                    gpvs: dict = read_gradient_property_vector_rpa(wf = self._gp, property = property,
                                        verbose = self._verbose, multiplicity=multiplicity[count])
                    if self._allr == None:
                        self._allr = False

            #Principal Propagator
            # - two integrals
            if not self._tintegrals:
                coulomb_integrals, exchange_integrals = get_coulomb_exchange_integrals(self._wf,
                                                        time_object = driver_time,
                                                        verbose = self._verbose,
                                                        verbose_int = verbose_integrals)
            else:
                raise ValueError("Falta implementar que obtenga las integrales desde 2i dadas")
            # - Build PP
            if response_type[str(len(response))] == "lineal":
                if not self._moe:
                    self._moe = self._wf.mo_energies
                principal_propagator: np.array = get_principal_propagator_lineal_rpa(
                                                                n_mo_occ = self._wf.mo_occ, n_mo_virt = self._wf.mo_virt,
                                                                moe = np.array(self._moe), coulomb = coulomb_integrals,
                                                                exchange = exchange_integrals,
                                                                multiplicity = multiplicity[0], tp_inv = 0,
                                                                time_object = driver_time,
                                                                verbose = self._verbose)
            else:
                raise ValueError("***ERROR \n\n\
                        Only is implemeted principal propagator of lineal response\n\
                        at rpa level.")

            #Response calculate
            calculate_lineal_reponse(n_mo_occ = self._wf.mo_occ, n_mo_virt = self._wf.mo_virt,
            principal_propagator = principal_propagator, gpvs = gpvs, all_responses = self._allr,
            time_object = driver_time, verbose = self._verbose)

            if self._verbose > 10:
                driver_time.add_name_delta_time(name = "Response Calculation", delta_time = (time() - start))
                driver_time.printing()
            print_title(name = f"END REPONSE CALCULATION")

if __name__ == "__main__":
    wfn = wave_function("../tests/molden_file/LiH.molden")
    r = response(wfn, properties = [["fc","fc"]], multiplicity=[[3,3]], verbose=12)
    r.rpa()