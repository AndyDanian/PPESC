from shielding import *
from libp import *

class ppesc():
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
    def drv_ppesc(self, principal_propagator_approximation: str = "rpa",
                    atoms: list = None,
                    ppesc_amounts: list = None,
                    isotropic: bool = True,
                    verbose: int = 0, verbose_integrals: int = -1):
        """
        Manage the PPESC calculation
        """

        print_title(name = "PPESC")

        if verbose > 10:
            start = time()
            driver_time = drv_time()
        else:
            driver_time = None

        if atoms is None:
            atoms = [*range(self._wf.atom_number)]

        # NR and Corrections
        all_responses, all_averages = \
            run_shielding_lresc(wf = self._wf, ppesc_amounts = ppesc_amounts,
                                    atom = atoms, ppesc_consts = ppesc_constants,
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
            driver_time.add_name_delta_time(name = "PPESC", delta_time = (time() - start))
            driver_time.printing()
        print_title(name = f"END LRESC CALCULATION")

if __name__ == "__main__":
    wfn = wave_function("../../tests/molden_file/HCl_v2z.molden")
    lr = ppesc(wfn)
    lr.drv_ppesc(verbose=11) #, ppesc_amounts = ["dianr"])