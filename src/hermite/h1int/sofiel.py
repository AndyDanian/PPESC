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
    if sofiel_xx and not (integrals._hermite_ao1b_binary.exists() and
                          integrals.binary(file = integrals._hermite_ao1b_binary, 
                                          label = "sofiel xx",
                                          io = "f")):
        symmetries["sofiel xx"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += charge[a] * integrals.binary(file=integrals._hermite_ao1b_binary,
                                                label=("nstcgo " + str(1 + a*3) + " x"), io="r")
        integrals.binary(file=integrals._hermite_ao1b_binary,dictionary={"sofiel xx": matrix},io="a")
        if verbose > 10:
            driver_time.add_name_delta_time(name = "SOFIEL XX AO", delta_time = (time() - start))

    if sofiel_xy and not (integrals._hermite_ao1b_binary.exists() and
                          integrals.binary(file = integrals._hermite_ao1b_binary, 
                                          label = "sofiel xy",
                                          io = "f")):
        symmetries["sofiel xy"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += charge[a] * integrals.binary(file=integrals._hermite_ao1b_binary,
                                                label=("nstcgo " + str(1 + a*3) + " y"), io="r")
        integrals.binary(file=integrals._hermite_ao1b_binary,dictionary={"sofiel xy": matrix},io="a")
        if verbose > 10:
            driver_time.add_name_delta_time(name = "SOFIEL XY AO", delta_time = (time() - start))

    if sofiel_xz and not (integrals._hermite_ao1b_binary.exists() and
                          integrals.binary(file = integrals._hermite_ao1b_binary, 
                                          label = "sofiel xz",
                                          io = "f")):
        symmetries["sofiel xz"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += charge[a] * integrals.binary(file=integrals._hermite_ao1b_binary,
                                                label=("nstcgo " + str(1 + a*3) + " z"), io="r")
        integrals.binary(file=integrals._hermite_ao1b_binary,dictionary={"sofiel xz": matrix},io="a")
        if verbose > 10:
            driver_time.add_name_delta_time(name = "SOFIEL XZ AO", delta_time = (time() - start))

    if sofiel_yx and not (integrals._hermite_ao1b_binary.exists() and
                          integrals.binary(file = integrals._hermite_ao1b_binary, 
                                          label = "sofiel yx",
                                          io = "f")):
        symmetries["sofiel yx"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += charge[a] * integrals.binary(file=integrals._hermite_ao1b_binary,
                                                label=("nstcgo " + str(2 + a*3) + " x"), io="r")
        integrals.binary(file=integrals._hermite_ao1b_binary,dictionary={"sofiel yx": matrix},io="a")
        if verbose > 10:
            driver_time.add_name_delta_time(name = "SOFIEL YX AO", delta_time = (time() - start))

    if sofiel_yy and not (integrals._hermite_ao1b_binary.exists() and
                          integrals.binary(file = integrals._hermite_ao1b_binary, 
                                          label = "sofiel yy",
                                          io = "f")):
        symmetries["sofiel yy"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += charge[a] * integrals.binary(file=integrals._hermite_ao1b_binary,
                                                label=("nstcgo " + str(2 + a*3) + " y"), io="r")
        integrals.binary(file=integrals._hermite_ao1b_binary,dictionary={"sofiel yy": matrix},io="a")
        if verbose > 10:
            driver_time.add_name_delta_time(name = "SOFIEL YY AO", delta_time = (time() - start))

    if sofiel_yz and not (integrals._hermite_ao1b_binary.exists() and
                          integrals.binary(file = integrals._hermite_ao1b_binary, 
                                          label = "sofiel yz",
                                          io = "f")):
        symmetries["sofiel yz"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += charge[a] * integrals.binary(file=integrals._hermite_ao1b_binary,
                                                label=("nstcgo " + str(2 + a*3) + " z"), io="r")
        integrals.binary(file=integrals._hermite_ao1b_binary,dictionary={"sofiel yz": matrix},io="a")
        if verbose > 10:
            driver_time.add_name_delta_time(name = "SOFIEL YZ AO", delta_time = (time() - start))

    if sofiel_zx and not (integrals._hermite_ao1b_binary.exists() and
                          integrals.binary(file = integrals._hermite_ao1b_binary, 
                                          label = "sofiel zx",
                                          io = "f")):
        symmetries["sofiel zx"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += charge[a] * integrals.binary(file=integrals._hermite_ao1b_binary,
                                                label=("nstcgo " + str(3 + a*3) + " x"), io="r")
        integrals.binary(file=integrals._hermite_ao1b_binary,dictionary={"sofiel zx": matrix},io="a")
        if verbose > 10:
            driver_time.add_name_delta_time(name = "SOFIEL ZX AO", delta_time = (time() - start))

    if sofiel_zy and not (integrals._hermite_ao1b_binary.exists() and
                          integrals.binary(file = integrals._hermite_ao1b_binary, 
                                          label = "sofiel zy",
                                          io = "f")):
        symmetries["sofiel zy"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += charge[a] * integrals.binary(file=integrals._hermite_ao1b_binary,
                                                label=("nstcgo " + str(3 + a*3) + " y"), io="r")
        integrals.binary(file=integrals._hermite_ao1b_binary,dictionary={"sofiel zy": matrix},io="a")
        if verbose > 10:
            driver_time.add_name_delta_time(name = "SOFIEL ZY AO", delta_time = (time() - start))

    if sofiel_zz and not (integrals._hermite_ao1b_binary.exists() and
                          integrals.binary(file = integrals._hermite_ao1b_binary, 
                                          label = "sofiel zz",
                                          io = "f")):
        symmetries["sofiel zz"] = "antisym"
        matrix = 0.0
        for a in range(number_atoms):
            matrix += charge[a] * integrals.binary(file=integrals._hermite_ao1b_binary,
                                                    label=("nstcgo " + str(3 + a*3) + " z"), io="r")
        integrals.binary(file=integrals._hermite_ao1b_binary,dictionary={"sofiel zz": matrix},io="a")
        if verbose > 10:
            driver_time.add_name_delta_time(name = "SOFIEL ZZ AO", delta_time = (time() - start))
