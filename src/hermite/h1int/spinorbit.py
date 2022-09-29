from lib1h import *

def spin_orbit(integrals: dict = None, number_atoms: int = None, charge: list = None,
                nprim: int = None,
                spino_x: bool = False, spino_y: bool = False, spino_z: bool = False,
                driver_time: object = None, verbose: int = 0):
    """
    Get spin orbit from PSO AO integrals

    Args:
        integrals (dict): Dictionary with all integrals calculated
        number_atoms (int): Number of atom
        charge (list): Atomic charge
        nprim (int): Number of primitives
        driver_time (object): Object to manage the time
        verbose (int): Print level
        spino_x (bool): Activate the calculate of spin orbit x-component
        spino_y (bool): Activate the calculate of spin orbit y-component
        spino_z (bool): Activate the calculate of spin orbit z-component

    """
    start = time()

    spinorbit_integrals: dict = {}
    symmetries: dict = {}
    matrix = np.zeros((nprim,nprim),dtype=float)

    SPINORBIT_CTE = 1.0/(137.0359998*137.0359998*4.0359998)

    if spino_x:
        symmetries["spinorbit x"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += SPINORBIT_CTE * charge[a] * integrals.binary(file=integrals._hermite_ao1b_binary,label=("pso " + str(1 + a*3)), io="r")
        integrals.binary(file=integrals._hermite_ao1b_binary,dictionary={"spinorbit x": matrix},io="a")
        spinorbit_integrals["spinorbit x"] = matrix
        if verbose > 10:
            driver_time.add_name_delta_time(name = "Spin-Orbit X AO", delta_time = (time() - start))

    if spino_y:
        symmetries["spinorbit y"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += SPINORBIT_CTE * charge[a] * integrals.binary(file=integrals._hermite_ao1b_binary,label=("pso " + str(2 + a*3)), io="r")
        spinorbit_integrals["spinorbit y"] = matrix
        if verbose > 10:
            driver_time.add_name_delta_time(name = "Spin-Orbit Y AO", delta_time = (time() - start))

    if spino_z:
        symmetries["spinorbit z"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += SPINORBIT_CTE * charge[a] * integrals.binary(file=integrals._hermite_ao1b_binary,label=("pso " + str(3 + a*3)), io="r")
        spinorbit_integrals["spinorbit z"] = matrix
        if verbose > 10:
            driver_time.add_name_delta_time(name = "Spin-Orbit Z AO", delta_time = (time() - start))
