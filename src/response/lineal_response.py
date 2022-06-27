from libr import *

def calculate_lineal_reponse(n_mo_occ: int = None, n_mo_virt: int = None,
                            principal_propagator: np.array = None, gpvs: dict = None,
                            all_responses: bool = False, verbose: int = 0):
    """
    Calculate of the path and total value of lienal response

    Args:
    ----
    n_mo_occ (int): Ocuppied molecular orbitals
    n_mo_virt (int): Virtual molecular orbitals
    moe (np.array, 1d): Molecular orbital energies
    coulomb (np.arra, 4d): Coulomb integrals
    exchange (np.array): Exchange integrals
    multiplicity (str): Multiplicity response
    tp_inv (int): Type of inverse: 0/numpy or 1/series
    all_response (bool): Activate the different response among one set
                        of operators
                        Example, if is True:
                        [A,B] then is get <<A;A>>, <<A;B>>, and <<B,B>>
    verbose (int): Print level
    """
    rotations: int = n_mo_occ*n_mo_virt
    names_gpv: list = [name for name in gpvs.keys()]
    values_gpv: list = [v for v in gpvs.values()]
    n_gpv: int = len(names_gpv)

    start = time()

    for index_gpv_left in range(n_gpv):
        for index_gpv_right in range(index_gpv_left, n_gpv):
            if not all_responses and index_gpv_left == index_gpv_right:
                continue

            if verbose > 20:
                title(name = f"PATHS: <<{names_gpv[index_gpv_left]},{names_gpv[index_gpv_right]}>>")

            ipath: int = 0
            vpathT: float = 0.0E+0
            irow: int = 0
            icol: int = 0
            for i in range(n_mo_occ):
                for a in range(n_mo_virt):
                    s = a + n_mo_occ

                    count = 0
                    spath = 0.0E+0
                    for j in range(n_mo_occ):
                        for b in range(n_mo_virt):
                            t = b + n_mo_occ

                            appb =\
                                values_gpv[index_gpv_left][a+i*n_mo_virt]\
                                *values_gpv[index_gpv_right][b+j*n_mo_virt]
                            # A_{i,s} PP_{i,s,j,t} B_{t,j}
                            appb = appb*principal_propagator[irow,icol]

                            vpathT += appb

                            if verbose > 20 and count == 0:
                                print(f" # ".center(8),f"i".center(8),
                                    f"s".center(8),f"j".center(8),
                                    f"t".center(8))
                            if verbose > 20 and abs(appb) > 0.1:
                                print(f"{count + 1}".center(8),f"{i}".center(8),
                                    f"{s}".center(8),f"{j}".center(8),
                                    f"{t}".center(8),f"{appb:.6f}".center(16))

                            spath += appb
                            count += 1
                            ipath += 1

                            icol +=1
                            if icol == rotations: irow += 1
                            if icol == rotations: icol = 0
                    if verbose > 20:
                        print("-"*60)
                        print(f'Total {spath:.6f}')
                        print()

            print()
            print(('='*40).center(70))
            print(f'-<<{names_gpv[index_gpv_left]};{names_gpv[index_gpv_right]}>>  =  {-vpathT:.6f}'.center(70))
            print(('='*40).center(70))

    if verbose > 10:
        name = f"Lineal Response"
        print_time(name = name,  delta_time=(time() - start))