from paramagnetic import run_shielding_paramagnetic
from diamagnetic import *
from libl import *

class lresc():
    def __init__(self, wf: wave_function = None):
        """
        Object to calculate LRESC values of different properties

        Args:
        ----
        wf (wave_function): Wave Function object
        """

        if not wf:
            raise ValueError("** ERROR \n\n\
                    The information is not enough to calculate. It is necessary\n\
                    wave function object.\
                    ")

        self._wf = wf
    ################################################################################################
    # METHODS
    ################################################################################################
    def drv_lresc(self, principal_propagator_approximation: str = "rpa",
                    atom: list = None, lresc_constant: str = "lresc",
                    paramagnetic: list = None, diamagnetic: list = None,
                    verbose: int = 0, verbose_responses: int = -1):
        """
        Manage the LRESC calculation
        """

        print_title(name = "LRESC")

        lresc_consts = lresc_constants[lresc_constant]

        if verbose > 10:
            start = time()
            driver_time = drv_time()
        else:
            driver_time = None

        if atom is None:
            atom = [*range(self._wf.atom_number)]

        # PARA: NR and Corrections
        #para_values = 
        run_shielding_paramagnetic(wf = self._wf, paramagnetic = paramagnetic,
                                    atom = atom, lresc_consts = lresc_consts,
                                    principal_propagator_approximation = principal_propagator_approximation,
                                    driver_time = driver_time, verbose = verbose,
                                    verbose_responses = verbose_responses)
        # DIA: NR and Corrections
        # dia_values = run_shielding_diamagnetic(wf = self._wf, diamagnetic = diamagnetic,
        #                                     atom = atom,
        #                                     principal_propagator_approximation = principal_propagator_approximation,
        #                                     driver_time = driver_time, verbose = verbose)

        if verbose > 10:
            driver_time.add_name_delta_time(name = "LRESC", delta_time = (time() - start))
            driver_time.printing()
        print_title(name = f"END LRESC CALCULATION")

if __name__ == "__main__":
    wfn = wave_function("../../tests/molden_file/LiH_pople.molden")
    lr = lresc(wfn)
    lr.drv_lresc(verbose=11, verbose_responses = 11, lresc_constant = "lresc_scale")