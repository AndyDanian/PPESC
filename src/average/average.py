
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
        verbose_integrals (int): Print level for integrals calculation
        """

        driver_time = self._wf._driver_time
        start = time()

        print_title("Average Value")

        # atomic integrals
        calculate_integral = eint(self._wf)
        integrals_1b, symmetries_1b = calculate_integral.integration_onebody(
        integrals_names = property, gauge = gauge, verbose = verbose_integrals)

        # molecular integrals
        n_mo_occ = self._wf.mo_occ
        mo_coeff_T = np.array(self._wf.mo_coefficients) #from molden come mo coefficient transpose

        mo_integral: dict = {}
        averages: dict = {}
        for name, atomic_int in integrals_1b.items():
            mo_integral[name] = np.matmul(mo_coeff_T,np.matmul(np.array(atomic_int), mo_coeff_T.T))
            averages[name] = 2.0*sum([mo_integral[name][i][i] for i in range(n_mo_occ)])

        if verbose >= 0:
            for name, average in averages.items():
                print_result(name = f"Average {name.title()}",
                value = f"{average:.8f}")

        driver_time.add_name_delta_time(name = f"Average Value", delta_time = (time() - start))
        driver_time.printing()
        driver_time.reset

        print_title("End Average Value")

        return averages

if __name__ == "__main__":
    wfn = wave_function("../tests/molden_file/LiH.molden")
    av = average(wfn)
    av.calculate_average(property = ["fc", "massvelo", "darwin"], verbose = 11, verbose_integrals=21)
