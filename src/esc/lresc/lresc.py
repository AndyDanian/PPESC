from typing import Union

from shielding import *
from libl import *


class lresc:
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
    def drv_lresc(
        self,
        atoms: list = [],
        lresc_amounts: list = [],
        principal_propagator_approximation: str = "rpa",
        lresc_constant: str = "lresc",
        tensor: bool = False,
        ani_axe: Union[str, int] = "z",
        verbose: int = 1,
        verbose_response: int = -1,
        verbose_average: int = -1,
    ):
        """
        Manage the LRESC calculation
        """
        ## Instance external objects
        # - Scratch
        io = self._wf._driver_scratch
        # - Driver Time
        driver_time = self._wf._driver_time
        ##
        start: float = time()

        io.write_output("LRESC", type=1)

        lresc_consts: dict[str, float] = lresc_constants[lresc_constant]

        if len(atoms) == 0:
            atoms = [*range(self._wf.atom_number)]
        # NR and Corrections
        all_responses, all_averages = run_shielding_lresc(
            io=io,
            driver_time=driver_time,
            verbose=verbose,
            wf=self._wf,
            lresc_amounts=lresc_amounts,
            atom=atoms,
            lresc_consts=lresc_consts,
            principal_propagator_approximation=principal_propagator_approximation,
            tensor=tensor,
            verbose_response=verbose_response,
            verbose_average=verbose_average,
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
            io.write_output(f"Isotropic Values [ppm]".center(76))
            io.write_output(f"and".center(76))
            io.write_output(f"Anisotropic Values [ppm]".center(76))

        io.write_output(("-" * 40).center(76))
        io.write_output(f"Anisotropic is with respect {ani_axe} axe.")
        io.write_output(
            "Values were multiplied by respectively constants according LRESC theory."
        )
        io.write_output("Paramagnetic corrections:")

        io.write_output(" Lineal Response:")
        io.write_output(
            "     * Singlet: {}, {}".format(
                lresc_label["lpsokin"], lresc_label["lkinpso"]
            )
        )
        io.write_output(
            "     * Triplet: {}, {}, {}, {}".format(
                lresc_label["fclap"],
                lresc_label["sdlap"],
                lresc_label["fcbso"],
                lresc_label["sdbso"],
            )
        )

        io.write_output(" Quadratic Response:")
        io.write_output(
            "     * Singlet: {}, {} ".format(
                lresc_label["lpsomv"], lresc_label["lpsodw"]
            )
        )
        io.write_output(
            "     * Triplet: {}, {}".format(lresc_label["lfcso"], lresc_label["lsdso"])
        )

        io.write_output("Diamagnetic corrections:")
        io.write_output(
            "      * Lineal Response: {}, {} ".format(
                lresc_label["a2mv"], lresc_label["a2dw"]
            )
        )

        io.write_output(
            "      *Averages: {}, {}, {} ".format(
                lresc_label["fc"], lresc_label["psooz"], lresc_label["dnske"]
            )
        )
        # Print results
        if not tensor:
            print_lresc_values(
                io=io,
                isotropic_responses=isotropic_responses,
                isotropic_averages=isotropic_averages,
                anisotropic_responses=anisotropic_responses,
                anisotropic_averages=anisotropic_averages,
                atom_label=self._wf.atomic_symbols,
                verbose=verbose,
            )
        else:
            print_lresc_tensor(
                io=io,
                isotropic_responses=isotropic_responses,
                isotropic_averages=isotropic_averages,
                anisotropic_responses=anisotropic_responses,
                anisotropic_averages=anisotropic_averages,
                responses_tensor=all_responses,
                averages_tensor=all_averages,
                atom_label=self._wf.atomic_symbols,
            )

        driver_time.add_name_delta_time(name="LRESC", delta_time=(time() - start))
        io.write_output(type=2, drv_time=driver_time)
        driver_time.reset
        io.write_output(f"END LRESC CALCULATION", type=1)


if __name__ == "__main__":
    wfn = wave_function(
        "../../tests/molden_file/H2_STO2G.molden",
        scratch_path="/home1/scratch",
        job_folder="160922134451",
    )
    lr = lresc(wfn)
    lr.drv_lresc(
        verbose=11,
        lresc_constant="lresc_scale",
        tensor=True,
        verbose_average=1,
        verbose_response=1,
    )
