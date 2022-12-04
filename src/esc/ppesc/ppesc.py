from typing import Union

from shielding import *
from libp import *


class ppesc:
    def __init__(self, wf: wave_function):
        """
        Object to calculate LRESC values of different properties

        Args:
        ----
        wf (wave_function): Wave Function object
        """

        if not wf:
            raise ValueError(
                "** ERROR \n\n\
                    The information is not enough to calculate. It is necessary\n\
                    wave function object.\
                    "
            )

        self._wf: wave_function = wf

    ################################################################################################
    # METHODS
    ################################################################################################
    def drv_ppesc(
        self,
        principal_propagator_approximation: str = "rpa",
        atoms: list = [],
        ppesc_amounts: list = [],
        tensor: bool = False,
        scalar_correction: bool = True,
        ani_axe: Union[str, int] = "z",
        verbose: int = 1,
        verbose_response: int = -1,
        verbose_average: int = -1,
        verbose_fock: int = 1,
    ):
        """
        Manage the PPESC calculation
        """
        ## Instance external objects
        # - Scratch
        io = self._wf._driver_scratch
        # - Driver Time
        driver_time = self._wf._driver_time
        ##
        start: float = time()
        io.write_output(information="PPESC", type=1)

        if len(atoms) == 0:
            atoms = [*range(self._wf.atom_number)]

        # NR and Corrections
        all_responses, all_averages = run_shielding(
            io=io,
            driver_time=driver_time,
            verbose=verbose,
            wf=self._wf,
            atom=atoms,
            ppesc_amounts=ppesc_amounts,
            ppesc_consts=ppesc_constants,
            scalar_correction=scalar_correction,
            principal_propagator_approximation=principal_propagator_approximation,
            tensor=tensor,
            verbose_response=verbose_response,
            verbose_average=verbose_average,
            verbose_fock=verbose_fock,
        )

        # Isotropoic and Anisitropic Value
        (
            isotropic_responses,
            isotropic_averages,
            anisotropic_responses,
            anisotropic_averages,
        ) = get_shielding_iso_ani(
            all_responses=all_responses,
            all_averages=all_averages,
            tensor=tensor,
            ani_axe=ani_axe,
        )
        if tensor:
            io.write_output(f"Tensor Values [ppm]".center(101))
        else:
            io.write_output(f"Isotropic Values [ppm]".center(101))
            io.write_output("and".center(101))
            io.write_output("Anisotropic Values [ppm]".center(101))
        io.write_output(("-" * 40).center(101))
        io.write_output(f"Anisotropic is with respect {ani_axe} axe.")
        io.write_output(
            "Values were multiplied by respectively constants according PPESC theory."
        )
        io.write_output("Paramagnetic corrections:")
        io.write_output(
            "     * Singlet: {}, {}".format(
                ppesc_label["lpsokin"], ppesc_label["lkinpso"]
            )
        )
        io.write_output(
            "     * Triplet: {}, {}".format(ppesc_label["fclap"], ppesc_label["sdlap"])
        )
        io.write_output("Diamagnetic corrections:")
        io.write_output(
            "      *Averages: {}, {}, {}, {}, {}".format(
                ppesc_label["fc"],
                ppesc_label["sd"],
                ppesc_label["psooz"],
                ppesc_label["dnske"],
                ppesc_label["pnstcgop"],
            )
        )

        if not tensor:
            # Print results
            print_ppesc_values(
                io=io,
                isotropic_responses=isotropic_responses,
                isotropic_averages=isotropic_averages,
                anisotropic_responses=anisotropic_responses,
                anisotropic_averages=anisotropic_averages,
                atom_label=self._wf.atomic_symbols,
                verbose=verbose,
            )
        else:
            print_ppesc_tensor(
                io=io,
                responses_tensor=all_responses,
                averages_tensor=all_averages,
                isotropic_responses=isotropic_responses,
                isotropic_averages=isotropic_averages,
                anisotropic_responses=anisotropic_responses,
                anisotropic_averages=anisotropic_averages,
                atom_label=self._wf.atomic_symbols,
            )

        driver_time.add_name_delta_time(name="PPESC", delta_time=(time() - start))
        io.write_output(type=2, drv_time=driver_time)
        driver_time.reset
        io.write_output(information=f"END LRESC CALCULATION", type=1)


if __name__ == "__main__":
    wfn = wave_function(
        "../../tests/molden_file/HI_v2z.molden",
        scratch_path="/home1/scratch",
        job_folder="HFSRv2z",
    )
    lr = ppesc(wfn)
    lr.drv_ppesc(
        verbose=11,
        scalar_correction=True,
        tensor=True,
        verbose_response=0,
        verbose_average=0,
        verbose_fock=11
    )
