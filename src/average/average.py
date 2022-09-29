
from liba import *

class average():
    def __init__(self, wf: wave_function = None):
        """
        Driver to calculate averge value

        Args:
            wf (wave_function): Wave function object
        """

        self._wf = wf

        if not wf:
            raise ValueError("** ERROR \n\n\
                    The information is not enogh to calculate. It is necessary\n\
                    wave function object.\
                    ")
    ################################################################################################
    # METHODS
    ################################################################################################
    def calculate_average(self, property: list = None, gauge: list = None, verbose: int = 0, verbose_integrals: int = -1):
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
        start = time()

        io.write_output(information = "Average Value", type = 1)

        # atomic integrals
        calculate_integral = eint(self._wf)
        calculate_integral.integration_onebody(
        integrals_names = property, gauge = gauge, verbose = verbose_integrals)

        # molecular integrals
        n_mo_occ = self._wf.mo_occ
        mo_coeff_T = np.array(self._wf.mo_coefficients) #from molden come mo coefficient transpose

        averages: dict = {}
        for label in calculate_integral.list_1b_integrals:
            mo_integral: np.array = np.matmul(  mo_coeff_T,
                                                np.matmul(
                                                io.binary(file = io._hermite_ao1b_binary, 
                                                          io = "r",
                                                          label = label),
                                                mo_coeff_T.T))
            
            if verbose >= 10:
                io.write_output(type = 9, direct = True, dictionary = {label: mo_integral})
            averages[label] = 2.0*sum([mo_integral[i,i] for i in range(n_mo_occ)])


        for name, average in averages.items():
            io.write_output(information = f"<{name.title()}>: {average:.8f}", type = 1, title_type = 2)

        driver_time.add_name_delta_time(name = f"Average Value", delta_time = (time() - start))
        io.write_output(drv_time = driver_time, type = 2)
        driver_time.reset

        io.write_output(information = "End Average Value", type = 1)

        return averages

if __name__ == "__main__":
    wfn = wave_function("../tests/molden_file/LiH.molden", scratch_path = "/home1/scratch", job_folder = "160922134451")
    av = average(wfn)
    av.calculate_average(property = ["fc", "massvelo", "darwin"], verbose = 11)
