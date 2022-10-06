from lib1h import *


def sofiel(
    integrals: scratch,
    driver_time: drv_time,
    number_atoms: int,
    charge: list[float],
    nprim: int,
    list_sofiel: list[bool],
    verbose: int = 0,
):
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
    start: float = time()
    sofiel_name: list[str] = [
        "sofiel xx",
        "sofiel xy",
        "sofiel xz",
        "sofiel yx",
        "sofiel yy",
        "sofiel yz",
        "sofiel zx",
        "sofiel zy",
        "sofiel zz",
    ]
    nstcgo_component = [
        lambda a: str(1 + a * 3) + " x",
        lambda a: str(1 + a * 3) + " y",
        lambda a: str(1 + a * 3) + " z",
        lambda a: str(2 + a * 3) + " x",
        lambda a: str(2 + a * 3) + " y",
        lambda a: str(2 + a * 3) + " z",
        lambda a: str(3 + a * 3) + " x",
        lambda a: str(3 + a * 3) + " y",
        lambda a: str(3 + a * 3) + " z",
    ]
    # * Calculation SOFIEL Integrals ****************************************
    for count, sofiel in enumerate(list_sofiel):
        if sofiel and not (
            integrals._hermite_ao1b_binary.exists()
            and integrals.binary(
                file=integrals._hermite_ao1b_binary, label=sofiel_name[count], io="f"
            )
        ):
            matrix: np.ndarray = np.zeros((nprim, nprim), dtype=float)
            for a in range(number_atoms):
                matrix += charge[a] * integrals.binary(
                    file=integrals._hermite_ao1b_binary,
                    label=("nstcgo " + nstcgo_component[count](a)),
                    io="r",
                )
            # Write: SOFIEL integrals in AO1BINT #################################
            integrals.binary(
                file=integrals._hermite_ao1b_binary,
                dictionary={sofiel_name[count]: matrix},
                io="a",
            )
            # END SOFIEL Integrals ################################################
            if verbose > 10:
                driver_time.add_name_delta_time(
                    name="SOFIEL " + nstcgo_component[count](a).upper() + " AO",
                    delta_time=(time() - start),
                )
    # * END Calculation of SOFIEL ******************************************************
