from libl import *

def correction_to_calculate(paramagnetic: list = None):
    """
    Turn on or off paramagnetic corrections
    Args:
        paramagneitc (list): Paramagnetic amount to calculate
    """
    if not paramagnetic:
        return {
        # Lineal -- Triplet
        "fckin": True, "sdkin": True, "fcbso": True, "sdbso": True,
        # Lineal -- Singlet
        "lpsokin": True, "lkinpso": True,
        # Quadratic -- Triplet
        "lsdspo": True,
        # Quadratic -- Singlet
        "lpsodw": True, "lpsomv": True
        }
    else:
        return {
        # Lineal -- Triplet
        "fckin": "fckin" in paramagnetic, "sdkin": "sdkin" in paramagnetic,
        "fcbso": "fcbso" in paramagnetic, "sdbso": "sdbso" in paramagnetic,
        # Lineal -- Singlet
        "lpsokin": "lpsokin" in paramagnetic, "lkinpso": "lkinpso" in paramagnetic,
        # Quadratic -- Triplet
        "lsdspo": "lsdspo" in paramagnetic,
        # Quadratic -- Singlet
        "lpsodw": "lpsodw" in paramagnetic, "lpsomv": "lpsomv" in paramagnetic
        }

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
        rps = response(wf = wf)
        gauge = wf.coordinates[a]
        corrections: dict = correction_to_calculate(paramagnetic)
        # -- Lineal Response
        # - Triplet
        # FcKin
        if corrections["fckin"]:
            fckin = rps.drv_reponse_calculation(principal_propagator_approximation = principal_propagator_approximation,
                                                properties = ["fc " + str(1 + a), "kinetic"], verbose = verbose)


    if verbose > 10:
        driver_time.add_name_delta_time(name = "Paramagnetic Calculations", delta_time = (time() - start))