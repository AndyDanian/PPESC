from io import BytesIO

from lib1h import *


def ozso(
    driver_time: drv_time,
    integrals: scratch,
    number_atoms: int,
    charge: list[float],
    nprim: int,
    list_spino: list[bool],
    verbose: int = 0,
):
    """
    Get spin orbit correction to orbital--zeeman from PSO-OZ AO integrals
    in spherical symmetry

    Args:
        integrals (dict): Dictionary with all integrals calculated
        number_atoms (int): Number of atom
        charge (list): Atomic charge
        nprim (int): Number of primitives
        driver_time (object): Object to manage the time
        verbose (int): Print level
        list_spino (bool): Activate the calculate of spin orbit:
                                                             x-component
                                                             y-component
                                                             z-component
    """
    start: float = time()
    ozso_name: list[str] = [
        "ozso x",
        "ozso y",
        "ozso z",
    ]
    psooz_component = [
        lambda a: str(1 + a * 3) + " x",
        lambda a: str(2 + a * 3) + " y",
        lambda a: str(3 + a * 3) + " z",
    ]
    #! SPINORBIT_CTE: float = 1.0 / (137.0359998 * 137.0359998 * 4.0359998) # DALTON
    SPINORBIT_CTE: float = 1.0 / (137.0359998 * 137.0359998 * 4.0)
    # * Calculation Spin-Orbit Correction to Orbital--Zeeman ***************************************
    for count, ozso in enumerate(list_spino):
        if ozso and not (
            integrals._hermite_ao1b_binary.exists()
            and integrals.binary(
                file=integrals._hermite_ao1b_binary, label=ozso_name[count], io="f"
            )
        ):
            matrix: np.ndarray = np.zeros((nprim, nprim), dtype=float)
            for a in range(number_atoms):
                matrix += (
                    SPINORBIT_CTE
                    * charge[a]
                    * integrals.binary(
                        file=integrals._hermite_ao1b_binary,
                        label=("psooz " + psooz_component[count](a)),
                        io="r",
                    )
                )
            # Write: Spin-Orbit Correction in AO1BINT #################################
            integrals.binary(
                file=integrals._hermite_ao1b_binary,
                dictionary={ozso_name[count]: matrix},
                io="a",
            )
            # END Spin-Orbit Correction ################################################
            if verbose > 10:
                driver_time.add_name_delta_time(
                    name="Spin-Orbit Correction to OZ", delta_time=(time() - start)
                )
    # * END Calculation of Spin-Orbit Correction to Orbital--Zeeman*********************************
