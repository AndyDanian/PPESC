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
                    tensor: bool = False,
                    scalar_correction: bool = True,
                    ani_axe: str or int = None,
                    verbose: int = 0, verbose_response: int = -1,
                    verbose_average: int = -1, verbose_fock: int = 1):
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
            run_shielding(wf = self._wf, ppesc_amounts = ppesc_amounts,
                            atom = atoms, ppesc_consts = ppesc_constants,
                            scalar_correction = scalar_correction,
                            principal_propagator_approximation = principal_propagator_approximation,
                            tensor = tensor,
                            driver_time = driver_time, verbose = verbose,
                            verbose_response = verbose_response,
                            verbose_average=verbose_average,
                            verbose_fock=verbose_fock)

        # Isotropoic and Anisitropic Value
        isotropic_responses, isotropic_averages,\
        anisotropic_responses, anisotropic_averages = (
                                        get_shielding_iso_ani(all_responses = all_responses,
                                        all_averages = all_averages, tensor = tensor,
                                        ani_axe = ani_axe)
                                        )
        if not tensor:
            # Print results
            print_ppesc_values(isotropic_responses = isotropic_responses,
                                isotropic_averages = isotropic_averages,
                                anisotropic_responses = anisotropic_responses,
                                anisotropic_averages = anisotropic_averages,
                                ani_axe = ani_axe,
                                atom_label = self._wf.atomic_symbols, verbose = verbose)
        else:
            print_ppesc_tensor(responses_tensor = all_responses, averages_tensor = all_averages,
                                isotropic_responses = isotropic_responses,
                                isotropic_averages = isotropic_averages,
                                anisotropic_responses = anisotropic_responses,
                                anisotropic_averages = anisotropic_averages, ani_axe = ani_axe,
                                atom_label = self._wf.atomic_symbols)

        if verbose > 10:
            driver_time.add_name_delta_time(name = "PPESC", delta_time = (time() - start))
            driver_time.printing()
        print_title(name = f"END LRESC CALCULATION")

if __name__ == "__main__":
    wfn = wave_function("../../tests/molden_file/HF_v2z.molden")
    lr = ppesc(wfn)
    lr.drv_ppesc(verbose=11, scalar_correction=False, tensor=True)