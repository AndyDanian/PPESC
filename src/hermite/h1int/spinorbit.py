from lib1h import *

def spin_orbit(integrals: dict = None, number_atoms: int = None, charge: list = None,
                spino_x: bool = False, spino_y: bool = False, spino_z: bool = False,
                driver_time: object = None, verbose: int = 0):
    """
    Get spin orbit from PSO AO integrals

    Args:
        integrals (dict): Dictionary with all integrals calculated
        number_atoms (int): Number of atom
        charge (list): Atomic charge
        driver_time (object): Object to manage the time
        verbose (int): Print level
        spino_x (bool): Activate the calculate of spin orbit x-component
        spino_y (bool): Activate the calculate of spin orbit y-component
        spino_z (bool): Activate the calculate of spin orbit z-component

    """
    start = time()

    spinorbit_integrals: dict = {}
    symmetries: dict = {}
    if spino_x:
        symmetries["spinorbit x"] = "antisym"
        for a in range(number_atoms):
            if a == 0:
                spinorbit_integrals["spinorbit x"] = [
                                charge[a] * value for value in integrals["pso " + str(1 + a*3)]
                                ]
            else:
                spinorbit_integrals["spinorbit x"] = [
                            old + charge[a] * new
                            for old, new in
                                zip(spinorbit_integrals["spinorbit x"], integrals["pso " + str(1 + a*3)])]
        if verbose > 10:
            driver_time.add_name_delta_time(name = "Spin-Orbit X AO", delta_time = (time() - start))
    if spino_y:
        symmetries["spinorbit y"] = "antisym"
        for a in range(number_atoms):
            if a == 0:
                spinorbit_integrals["spinorbit y"] = [
                                charge[a] * value for value in integrals["pso " + str(2 + a*3)]
                                ]
            else:
                spinorbit_integrals["spinorbit y"] = [
                            old + charge[a] * new
                            for old, new in
                                zip(spinorbit_integrals["spinorbit y"], integrals["pso " + str(2 + a*3)])]
        if verbose > 10:
            driver_time.add_name_delta_time(name = "Spin-Orbit Y AO", delta_time = (time() - start))
    if spino_z:
        symmetries["spinorbit z"] = "antisym"
        for a in range(number_atoms):
            if a == 0:
                spinorbit_integrals["spinorbit z"] = [
                                charge[a] * value for value in integrals["pso " + str(3 + a*3)]
                                ]
            else:
                spinorbit_integrals["spinorbit z"] = [
                            old + charge[a] * new
                            for old, new in
                                zip(spinorbit_integrals["spinorbit z"], integrals["pso " + str(3 + a*3)])]
        if verbose > 10:
            driver_time.add_name_delta_time(name = "Spin-Orbit Z AO", delta_time = (time() - start))

    return spinorbit_integrals, symmetries