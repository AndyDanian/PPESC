from liba import *


class average:
    def __init__(self, wf: wave_function):
        """
        Driver to calculate averge value

        Args:
            wf (wave_function): Wave function object
        """

        if not wf:
            raise ValueError(
                "** ERROR \n\n\
                    The information is not enogh to calculate. It is necessary\n\
                    wave function object.\
                    "
            )
        self._wf: wave_function = wf

    ################################################################################################
    # METHODS
    ################################################################################################
    def calculate_average(
        self,
        property: list[str] = [],
        gauge: list[float] = [],
        verbose: int = 0,
        verbose_integrals: int = -1,
    ) -> dict[str, float]:
        """
        Calculate reponse at random phase approximation

        Args:
        ----
        property (list): Properties name to calculate its average value
        gauge (list): Gauge origen
        verbose (int): Print level
                        * < 10: Minimum
                        * > 10: MO integrals
        verbose_integrals (int): Print level for integrals calculation
        """

        ## Instance external objects
        # - Scratch
        io = self._wf._driver_scratch
        # - Driver Time
        driver_time = self._wf._driver_time
        ##
        start: float = time()

        io.write_output(information="Average Value", type=1)

        # atomic integrals
        calculate_integral = eint(self._wf)
        calculate_integral.integration_onebody(
            integrals_names=property, gauge=gauge, verbose=verbose_integrals
        )

        # molecular integrals
        n_mo_occ: int = self._wf.mo_occ
        mo_coeff_T: np.ndarray = np.array(
            self._wf.mo_coefficients
        )  #  molden, mo coefficient matrix is transpose

        averages: dict[str, float] = {}
        for label in calculate_integral.list_1b_integrals:
            mo_integral: np.ndarray = np.matmul(
                mo_coeff_T,
                np.matmul(
                    io.binary(file=io._hermite_ao1b_binary, io="r", label=label),
                    mo_coeff_T.T,
                ),
            )

            if verbose >= 20:
                io.write_output(type=9, direct=True, dictionary={label: mo_integral})
            averages[label] = 2.0 * sum([mo_integral[i, i] for i in range(n_mo_occ)])

        for name, average in averages.items():
            io.write_output(
                information=f"<{name.title()}>: {average:.8f}", type=1, title_type=2
            )

        driver_time.add_name_delta_time(
            name=f"Average Value", delta_time=(time() - start)
        )
        io.write_output(drv_time=driver_time, type=2)
        driver_time.reset

        io.write_output(information="End Average Value", type=1)

        return averages


if __name__ == "__main__":
    wfn = wave_function(
        "../tests/molden_file/HAt_v2z.molden",
        scratch_path="/home1/scratch",
        #job_folder="HAtaverage",
        job_folder="/home1/build/PPESC/integrals/intpy/pySCF_CAL/HAt/DwMv/Av",
    )
    av = average(wfn)
    av.calculate_average(property=[
                                   "overlap",
                                   "nstcgo 1 x", "nstcgo 2 y", "nstcgo 3 z",
                                   "psooz 1 x", "psooz 2 y", "psooz 3 z",
                                   "fc 1", "sd 1 x", "sd 2 y", "sd 3 z",
                                   "dnske 1 x", "dnske 2 y", "dnske 3 z" 
                                   "cdnske 1 x", "cdnske 2 y", "cdnske 3 z"], 
                                   verbose=1)
