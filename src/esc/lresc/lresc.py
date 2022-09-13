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
                    tensor: bool = False,
                    isotropic: bool = True,
                    ani_axe: str or int = "z",
                    verbose: int = 1, verbose_response: int = -1,
                    verbose_average: int = -1):
        """
        Manage the LRESC calculation
        """

        print_title(name = "LRESC")

        lresc_consts = lresc_constants[lresc_constant]

        start = time()
        driver_time = drv_time()

        if atoms is None:
            atoms = [*range(self._wf.atom_number)]
        # NR and Corrections
        all_responses, all_averages = \
            run_shielding_lresc(wf = self._wf, lresc_amounts = lresc_amounts,
                                atom = atoms, lresc_consts = lresc_consts,
                                principal_propagator_approximation = principal_propagator_approximation,
                                driver_time = driver_time, verbose = verbose,
                                tensor = tensor,
                                verbose_response = verbose_response,
                                verbose_average = verbose_average)
        # Isotropoic and Anisitropic Value
        isotropic_responses, isotropic_averages,\
        anisotropic_responses, anisotropic_averages = get_shielding_iso_ani(
                                        all_responses = all_responses,
                                        all_averages = all_averages,
                                        tensor = tensor,
                                        ani_axe = ani_axe
                                        )
        if tensor:
            print(f"Tensor Values [ppm]".center(101))
        else:
            print(f"Isotropic Values [ppm]".center(76), "\n",
                    "and".center(76), "\n", "Anisotropic Values [ppm]".center(76))
        print(("-"*40).center(76))
        print("Anisotropic is with respect ",ani_axe," axe.")
        print("Values were multiplied by respectively constants according LRESC theory.")
        if verbose > 10 or tensor or len(atoms) < 3:
            print("Paramagnetic corrections:")
            print(" Lineal Response:")
            print("     * Singlet: ",lresc_label["lpsokin"],lresc_label["lkinpso"])
            print("     * Triplet: ",lresc_label["fclap"],lresc_label["sdlap"],
                                    lresc_label["fcbso"],lresc_label["sdbso"])
            print(" Quadratic Response:")
            print("     * Singlet: ",lresc_label["lpsomv"],lresc_label["lpsodw"])
            print("     * Triplet: ",lresc_label["lfcso"],lresc_label["lsdso"])
            print("Diamagnetic corrections:")
            print("      * Lineal Response: ",lresc_label["a2mv"],lresc_label["a2dw"])
            print("      *Averages: ",lresc_label["fc"],lresc_label["psooz"],
                    lresc_label["dnske"])
        # Print results
        if not tensor:
            print_lresc_values(isotropic_responses = isotropic_responses,
                                isotropic_averages = isotropic_averages,
                                anisotropic_responses = anisotropic_responses,
                                anisotropic_averages = anisotropic_averages,
                                atom_label = self._wf.atomic_symbols,
                                verbose = verbose)
        else:
            print_lresc_tensor(responses_tensor = all_responses, averages_tensor = all_averages,
                                isotropic_responses = isotropic_responses,
                                isotropic_averages = isotropic_averages,
                                anisotropic_responses = anisotropic_responses,
                                anisotropic_averages = anisotropic_averages,
                                atom_label = self._wf.atomic_symbols)

        if verbose > 0:
            driver_time.add_name_delta_time(name = "LRESC", delta_time = (time() - start))
            driver_time.printing()
        print_title(name = f"END LRESC CALCULATION")

if __name__ == "__main__":
    wfn = wave_function("../../tests/molden_file/LiH_pople.molden")
    lr = lresc(wfn)
    lr.drv_lresc(verbose=11, lresc_constant = "lresc_scale", tensor=False) #, lresc_amounts = ["lpsomv", "lpsodw", "lfcso", "lsdsoxx"])