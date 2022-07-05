from libl import *


def run_shielding_diamagnetic(wf: wave_function = None, diamagnetic: list = None,
                                atom: list = None,
                                principal_propagator_approximation: str = None,
                                driver_time: drv_time = None, verbose: int = 0):
    """
    Calculate all diamagnetic values
    Args:
        wf (wave_function): Wave function object
        diamagneitc (list): Diamagnetic amount to calculate
        atom (list): Atom index to calculate the shielding
        principal_propagator_approximation (str): Approximation
        driver_time (drv_time): Object to manage the time calculation
        verbose (int): Print level
    """
    start = time()

    if verbose > 10:
        driver_time.add_name_delta_time(name = "Diamagnetic Calculations", delta_time = (time() - start))
