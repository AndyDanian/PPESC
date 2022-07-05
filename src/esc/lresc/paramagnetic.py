from libl import *


def run_shielding_paramagnetic(wf: wave_function = None, paramagnetic: list = None,
                                atom: list = None,
                                principal_propagator_approximation: str = None,
                                driver_time: drv_time = None, verbose: int = 0):
    """
    Calculate all paramagnetic values
    Args:
        wf (wave_function): Wave function object
        paramagneitc (list): Paramagnetic amount to calculate
        atom (list): Atom index to calculate the shielding
        principal_propagator_approximation (str): Approximation
        driver_time (drv_time): Object to manage the time calculation
        verbose (int): Print level
    """
    start = time()

    for a in atom:
        gauge = wf.coordinates[a]
        fckin = ["fc " + str(1 + a) + ", kinenerg"]
        sd1kin = ["sd " + str(1 + a*3) + ", laplacian 1"] #In old LRESC used DPTOVL, which accoring of DALTON
        sd2kin = ["sd " + str(2 + a*3) + ", laplacian 2"] #source is -1/4*ddi
        sd3kin = ["sd " + str(3 + a*3) + ", laplacian 3"]


    if verbose > 10:
        driver_time.add_name_delta_time(name = "Paramagnetic Calculations", delta_time = (time() - start))