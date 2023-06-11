from liba import *


class projint:
    def __init__(self, wf: wave_function):
        """
        Driver to calculate projection integrals

        This module used different integrals to build one, used:
             <phi|A|Psi><Psi|B|phi>   phi: atomic orbital
                                      Psi: Molecular orbital
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
    def projections_integrals(
        self,
        integrals: list[list[str]] = [],
        gauge: list[float] = [],
        dipole: list[float] = [],
        verbose: int = 0,
        verbose_integrals: int = -1,
    ) -> dict[str, float]:
        """
        Calculate one integral from others

        Args:
        ----
        integrals (list): Name of integrals to build another. Set integrals 
                          for one are seperated by list
                          [[a, b],[c, d, e], ...]
        gauge (list):     Gauge origen
        r_dipole (list):  Dipole position
        verbose (int):    Print level
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

        for names in integrals:
            # Integrals will save in binary file
            calculate_integral.integration_onebody(
                integrals_names=names, gauge=gauge, dipole=dipole, verbose=verbose_integrals
            )

        #! molden, mo coefficient matrix is transpose
        mo_coeff_T: np.ndarray = np.array(
            self._wf.mo_coefficients
        )  
        # Build Projector
        projector: np.ndarray = np.matmul(mo_coeff_T.T,mo_coeff_T)
        if verbose >= 10:
            io.write_output(type=9, direct=True, dictionary={"|Psi><Psi|": projector})

        intprojectionint: dict[str, np.ndarray] = {}
        #for label in calculate_integral.list_1b_integrals:
        for names in integrals:
            label: str = ''
            for i, name in enumerate(names):
                if i == 0:
                    operator_a: np.ndarray() = (
                        io.binary(file=io._hermite_ao1b_binary, io="r", label=name)
                        #mo_coeff_T.T
                    )
                    label = name 
                else:
                    projection = np.matmul(
                                operator_a,
                                np.matmul(
                                projector,
                                io.binary(file=io._hermite_ao1b_binary, io="r", label=name)
                                ))
                    label += ', ' + name
            
            intprojectionint[label] = projection

            if verbose >= 10:
                io.write_output(type=9, direct=True, dictionary={label: projection})

        driver_time.add_name_delta_time(
            name=f"Projection on Integration <phi|A.1.B|phi>", delta_time=(time() - start)
        )
        io.write_output(drv_time=driver_time, type=2)
        driver_time.reset

        io.write_output(information="End <phi|A.1.B|phi>", type=1)

        return intprojectionint


if __name__ == "__main__":
    wfn = wave_function(
        "../tests/molden_file/HF_v2z.molden",
        scratch_path="/home1/scratch",
        job_folder="/home1/scratch/160922134451/",
    )
    p1p = projint(wfn)
    p1p.projections_integrals(integrals=[
                                         ["dxdx", "angmom x"],
                                         ["dydy", "angmom x"],
                                         ["dzdz", "angmom x"],
                                         ["angmom x", "dxdx"],
                                         ["angmom x", "dydy"],
                                         ["angmom x", "dzdz"],
                                         ["kinetic", "angmom x"],
                                         ["angmom x", "kinetic"],
                                        ],
                            verbose=11)
