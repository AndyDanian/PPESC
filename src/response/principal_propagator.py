from coulomb_exchange import *
from libr import *


def get_principal_propagator_lineal_rpa(
    io: scratch,
    driver_time: drv_time,
    n_mo_occ: int,
    n_mo_virt: int,
    moe: np.ndarray,
    multiplicity: str,
    tp_inv: int = 0,
    verbose: int = 0,
) -> np.ndarray:
    """
    Build principal propagator for lineal response at rpa level

    Args:
    ----
    io (object:scratch): Driver to driver the output and binary files
    n_mo_occ (int): Ocuppied molecular orbitals
    n_mo_virt (int): Virtual molecular orbitals
    moe (np.array, 1d): Molecular orbital energies
    coulomb (np.arra, 4d): Coulomb integrals
    exchange (np.array): Exchange integrals
    multiplicity (str): Multiplicity response
    tp_inv (int): Type of inverse: 0/numpy or 1/series
    quadratic (bool): Quadratic response
    verbose (int): Print level, 2D array

    Return:
    ------
    pp (np.array) : Principal propagator
    """
    start: float = time()

    D1: float = 1.0e0
    k: int = 0
    irow: int = 0
    icol: int = 0
    rotations: int = n_mo_virt * n_mo_occ

    if tp_inv == 1:
        fock_inv: np.ndarray = np.zeros((rotations, rotations), dtype=float)

    w: np.ndarray = np.zeros((rotations, rotations), dtype=float)
    coulomb: np.ndarray = io.binary(file=io._exchange_coulomb, io="r", label="coulomb")
    exchange: np.ndarray = io.binary(
        file=io._exchange_coulomb, io="r", label="exchange"
    )

    for i in range(n_mo_occ):
        for a in range(n_mo_virt):
            s: int = a + n_mo_occ
            if tp_inv == 1:
                fock_inv[k, k] = D1 / (moe[s] - moe[i])
            k += 1

            for j in range(n_mo_occ):
                for b in range(n_mo_virt):

                    delta_moe: float = 0.0
                    if i == j and a == b and tp_inv == 0:
                        delta_moe = moe[s] - moe[i]
                    if (
                        isinstance(multiplicity, str)
                        and multiplicity.lower() == "singlet"
                    ) or (isinstance(multiplicity, int) and multiplicity == 1):

                        w[irow, icol] = (
                            delta_moe + exchange[a, j, b, i] - coulomb[a, b, j, i]
                        )

                    elif (
                        isinstance(multiplicity, str)
                        and multiplicity.lower() == "triplet"
                    ) or (
                        isinstance(multiplicity, int) and multiplicity == 3
                    ):  # <ab|ji> + <aj|bi>

                        w[irow, icol] = (
                            delta_moe - coulomb[a, b, j, i] - exchange[a, j, b, i]
                        )

                    else:
                        raise ValueError(
                            "***ERROR\n\n\
                            Multiplicity can be 1/singlet or 3/triplet"
                        )

                    icol = icol + 1
                    if icol > rotations - 1:
                        irow += 1  # new row is equivalent to a change of i and a-index
                    if icol > rotations - 1:
                        icol = 0  #

    build_time: float = time() - start

    start = time()
    if tp_inv == 0:
        pp: np.ndarray = np.linalg.inv(w)
    else:
        raise ValueError(
            "***ERROR\n\n\
            Inverse of a matrix only with numpy is implemeted"
        )
    inverse_time: float = time() - start

    driver_time.add_name_delta_time(
        name=f"Build Principal Propagator", delta_time=build_time
    )
    driver_time.add_name_delta_time(name=f"Inverse", delta_time=inverse_time)

    if verbose > 50:
        label: str = f"Principal Propagator (Multiplicity: {multiplicity})"
        io.write_output(type=9, direct=True, dictionary={label: w})
        if tp_inv == 0:
            label: str = (
                f"Principal Propagator Inverse (Numpy, Multiplicity: {multiplicity})"
            )
        io.write_output(type=9, direct=True, dictionary={label: w})
    return pp


def drv_principal_propagator(
    io: scratch,
    driver_time: drv_time,
    moe: np.ndarray,
    n_mo_occ: int,
    n_mo_virt: int,
    multiplicity_pp: dict[str, bool],
    tp_inv: int,
    verbose: int = 0,
) -> None:
    """
    Driver to build the principal propagator

    Args:
    ----
    io (object:scratch): Driver to driver the output and binary files
    driver_time (drv_time): Object relationed with the manage of the time
    n_mo_occ (int): Ocuppied molecular orbitals
    n_mo_virt (int): Virtual molecular orbitals
    moe (list): Molecular orbital energies
    coulomb (np.array): Coulomb integrals in molecular base
    Exchange (np.array): Exchange integrals in molecular base
    multiplicity_pp (dict): Indicate if is neccesary calculate singlet or triplet principal
                            propagator
    tp_inv (int): Type of inverse: 0/numpy or 1/series
    verbose (int): Print level
    quadratic (bool): quadratic response
    verbose_integrals (int): Print level for hermite module
    """
    # - Build PP
    principal_propagator: dict[str, np.ndarray] = {}
    symmetries: str = ""
    for name, ms in multiplicity_pp.items():
        if ms:
            principal_propagator[name] = get_principal_propagator_lineal_rpa(
                io=io,
                driver_time=driver_time,
                n_mo_occ=n_mo_occ,
                n_mo_virt=n_mo_virt,
                moe=moe,
                multiplicity=name,
                tp_inv=tp_inv,
                verbose=verbose,
            )
            # Write in binary file the principal propagator
            io.binary(
                file=io._principal_propagator,
                dictionary={f"{name}": principal_propagator[name]},
                io="a",
            )
            symmetries += " " + name
    ## Write in the ouput the PRINPROP.H5 file
    io.write_output(
        information=f"Principal Propagator (symmetry: {symmetries})",
        type=3,
        size_file=io._principal_propagator.stat().st_size,
    )
