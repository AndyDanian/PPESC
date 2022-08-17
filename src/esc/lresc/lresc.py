from shielding import *
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
                    atoms: list = None, lresc_constant: str = "lresc",
                    lresc_amounts: list = None,
                    isotropic: bool = True,
                    verbose: int = 0, verbose_integrals: int = -1):
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

        if atoms is None:
            atoms = [*range(self._wf.atom_number)]

        # NR and Corrections
        all_responses, all_averages = \
            run_shielding_lresc(wf = self._wf, lresc_amounts = lresc_amounts,
                                    atom = atoms, lresc_consts = lresc_consts,
                                    principal_propagator_approximation = principal_propagator_approximation,
                                    driver_time = driver_time, verbose = verbose,
                                    verbose_integrals = verbose_integrals)

        # Isotropoic and Anisitropic Value
        isotropic_responses, isotropic_averages = get_shielding_isotropic(all_responses = all_responses,
                                        all_averages = all_averages)
        anisotropic_responses, anisotropic_averages = get_shielding_anisotropic(all_responses = all_responses,
                                        all_averages = all_averages)
        # Print results
        print_lresc_values(isotropic_responses = isotropic_responses, isotropic_averages = isotropic_averages,
                            anisotropic_responses = anisotropic_responses, anisotropic_averages = anisotropic_averages,
                            atom_label = self._wf.atomic_symbols)


        if verbose > 10:
            driver_time.add_name_delta_time(name = "LRESC", delta_time = (time() - start))
            driver_time.printing()
        print_title(name = f"END LRESC CALCULATION")

if __name__ == "__main__":
    wfn = wave_function("../../tests/molden_file/LiH_pople.molden")
    lr = lresc(wfn)
    lr.drv_lresc(verbose=1, lresc_constant = "lresc_scale") #, lresc_amounts = ["lpsomv", "lpsodw", "lfcso", "lsdsoxx"])