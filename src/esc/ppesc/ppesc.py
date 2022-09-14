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
                    ani_axe: str or int = "z",
                    verbose: int = 1, verbose_response: int = -1,
                    verbose_average: int = -1, verbose_fock: int = 1):
        """
        Manage the PPESC calculation
        """
        print_title(name = "PPESC")

        start = time()
        driver_time = drv_time()

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
        if tensor:
            print(f"Tensor Values [ppm]".center(101))
        else:
            print(f"Isotropic Values [ppm]".center(76), "\n", "and".center(76), "\n", "Anisotropic Values [ppm]".center(76))
        print(("-"*40).center(76))
        print("Anisotropic is with respect ",ani_axe," axe.")
        print("Values were multiplied by respectively constants according PPESC theory.")
        if tensor or len(atoms) < 3 or verbose > 10:
            print("Paramagnetic corrections:")
            print("     * Singlet: ",ppesc_label["lpsokin"],ppesc_label["lkinpso"])
            print("     * Triplet: ",ppesc_label["fclap"],ppesc_label["sdlap"])
            print("Diamagnetic corrections:")
            print("      *Averages: ",ppesc_label["fc"],ppesc_label["sd"],ppesc_label["psooz"],
                ppesc_label["dnske"], ppesc_label["pnstcgop"])

        if not tensor:
            # Print results
            print_ppesc_values(isotropic_responses = isotropic_responses,
                                isotropic_averages = isotropic_averages,
                                anisotropic_responses = anisotropic_responses,
                                anisotropic_averages = anisotropic_averages,
                                atom_label = self._wf.atomic_symbols, verbose = verbose)
        else:
            print_ppesc_tensor(responses_tensor = all_responses, averages_tensor = all_averages,
                                isotropic_responses = isotropic_responses,
                                isotropic_averages = isotropic_averages,
                                anisotropic_responses = anisotropic_responses,
                                anisotropic_averages = anisotropic_averages,
                                atom_label = self._wf.atomic_symbols)

        if verbose > 0:
            driver_time.add_name_delta_time(name = "PPESC", delta_time = (time() - start))
            driver_time.printing()
        print_title(name = f"END LRESC CALCULATION")

if __name__ == "__main__":
    wfn = wave_function("../../tests/molden_file/LiH_pople.molden")
    lr = ppesc(wfn)
    lr.drv_ppesc(verbose=11, scalar_correction=False, tensor=False, verbose_response=11, verbose_average=11)