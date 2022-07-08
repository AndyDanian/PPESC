from lib1h import *

def sofiel(integrals: dict = None, number_atoms: int = None, charge: list = None,
            nprim: int = None,
            sofiel_xx: bool = False, sofiel_yy: bool = False, sofiel_zz: bool = False,
            sofiel_xy: bool = False, sofiel_xz: bool = False, sofiel_yz: bool = False,
            sofiel_yx: bool = False, sofiel_zx: bool = False, sofiel_zy: bool = False,
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
        sofiel_xx (bool): Activate the calculate of spin orbit xx-component
        sofiel_yy (bool): Activate the calculate of spin orbit yy-component
        sofiel_zz (bool): Activate the calculate of spin orbit zz-component
        sofiel_xy (bool): Activate the calculate of spin orbit xy-component
        sofiel_xz (bool): Activate the calculate of spin orbit xz-component
        sofiel_yz (bool): Activate the calculate of spin orbit yz-component
        sofiel_yx (bool): Activate the calculate of spin orbit yx-component
        sofiel_zx (bool): Activate the calculate of spin orbit zx-component
        sofiel_zy (bool): Activate the calculate of spin orbit zy-component
    """
    start = time()

    spinorbit_integrals: dict = {}
    symmetries: dict = {}
    matrix = np.zeros((nprim,nprim),dtype=float)
    if sofiel_xx:
        symmetries["sofiel xx"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += charge[a] * integrals["nstcgo " + str(1 + a*3) + " x"]
        spinorbit_integrals["sofiel xx"] = matrix
        if verbose > 10:
            driver_time.add_name_delta_time(name = "SOFIEL XX AO", delta_time = (time() - start))

    if sofiel_xy:
        symmetries["sofiel xy"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += charge[a] * integrals["nstcgo " + str(1 + a*3) + " y"]
        spinorbit_integrals["sofiel xy"] = matrix
        if verbose > 10:
            driver_time.add_name_delta_time(name = "SOFIEL XY AO", delta_time = (time() - start))

    if sofiel_xz:
        symmetries["sofiel xz"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += charge[a] * integrals["nstcgo " + str(1 + a*3) + " z"]
        spinorbit_integrals["sofiel xz"] = matrix
        if verbose > 10:
            driver_time.add_name_delta_time(name = "SOFIEL XZ AO", delta_time = (time() - start))

    if sofiel_yx:
        symmetries["sofiel yx"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += charge[a] * integrals["nstcgo " + str(2 + a*3) + " x"]
        spinorbit_integrals["sofiel yx"] = matrix
        if verbose > 10:
            driver_time.add_name_delta_time(name = "SOFIEL YX AO", delta_time = (time() - start))

    if sofiel_yy:
        symmetries["sofiel yy"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += charge[a] * integrals["nstcgo " + str(2 + a*3) + " y"]
        spinorbit_integrals["sofiel yy"] = matrix
        if verbose > 10:
            driver_time.add_name_delta_time(name = "SOFIEL YY AO", delta_time = (time() - start))

    if sofiel_yz:
        symmetries["sofiel yz"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += charge[a] * integrals["nstcgo " + str(2 + a*3) + " z"]
        spinorbit_integrals["sofiel yz"] = matrix
        if verbose > 10:
            driver_time.add_name_delta_time(name = "SOFIEL YZ AO", delta_time = (time() - start))

    if sofiel_zx:
        symmetries["sofiel zx"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += charge[a] * integrals["nstcgo " + str(3 + a*3) + " x"]
        spinorbit_integrals["sofiel zx"] = matrix
        if verbose > 10:
            driver_time.add_name_delta_time(name = "SOFIEL ZX AO", delta_time = (time() - start))

    if sofiel_zy:
        symmetries["sofiel zy"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += charge[a] * integrals["nstcgo " + str(3 + a*3) + " y"]
        spinorbit_integrals["sofiel zy"] = matrix
        if verbose > 10:
            driver_time.add_name_delta_time(name = "SOFIEL ZY AO", delta_time = (time() - start))

    if sofiel_zz:
        symmetries["sofiel zz"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += charge[a] * integrals["nstcgo " + str(3 + a*3) + " z"]
        spinorbit_integrals["sofiel zz"] = matrix
        if verbose > 10:
            driver_time.add_name_delta_time(name = "SOFIEL ZZ AO", delta_time = (time() - start))


    return spinorbit_integrals, symmetries