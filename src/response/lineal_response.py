from libr import *

def calculate_lineal_reponse(operator_a: list = None, operator_b: list = None,
                            n_mo_occ: int = None, n_mo_virt: int = None,
                            principal_propagator: np.array = None, gpvs: dict = None,
                            time_object: drv_time = None,
                            verbose: int = 0):
    """
    Calculate of the path and total value of stactic lineal response

    Args:
    ----
    operator_a (list): Name of the operator in the left
    operator_b (list): Name of the operator in the right
    n_mo_occ (int): Ocuppied molecular orbitals
    n_mo_virt (int): Virtual molecular orbitals
    principal_propagator_a (np.array): Inverse of the principal propagator
    time_object (drv_time): Manage time calculation
    verbose (int): Print level
    """
    start = time()

    lineal_responses: dict = {}
    n_rotation: int = n_mo_occ*n_mo_virt
    for index_a, op_a in enumerate(operator_a):
        for index_b, op_b in enumerate(operator_b):
            if index_a > index_b and op_a.split()[0] == op_b.split()[0]:
                continue

            if verbose > 20:
                print_subtitle(name = f"PATHS: <<{op_a};{op_b}>>")

            ipath: int = 0
            vpathT: float = 0.0E+0
            for i in range(n_mo_occ):
                for a in range(n_mo_virt):
                    s = a + n_mo_occ

                    count = 0
                    spath = 0.0E+0
                    for j in range(n_mo_occ):
                        for b in range(n_mo_virt):
                            t = b + n_mo_occ

                            # -0.5*(<i|A|a>P_{ia,jb}^{-1}<b|B|j>
                            # +
                            # <i|B|b>P_{jb,ia}^{-1}<a|A|j>)

                            appb = -0.5*(
                                    gpvs[op_a][a+i*n_mo_virt]*
                                    principal_propagator[a+i*n_mo_virt,b+j*n_mo_virt]*
                                    gpvs[op_b][b+j*n_mo_virt+n_rotation]
                                +
                                    gpvs[op_b][b+j*n_mo_virt]*
                                    principal_propagator[b+j*n_mo_virt,a+i*n_mo_virt]*
                                    gpvs[op_a][a+i*n_mo_virt+n_rotation]
                                )

                            vpathT += appb

                            if verbose > 20 and count == 0:
                                print(f" # ".center(8),f"i".center(8),
                                    f"s".center(8),
                                    f"t".center(8),f"j".center(8))
                            if verbose > 20 and abs(appb) > 0.1:
                                print(f"{count + 1}".center(8),f"{i + 1}".center(8),
                                    f"{s + 1}".center(8),f"{t + 1}".center(8),
                                    f"{j + 1}".center(8),f"{appb:.6f}".center(16)
                                    )

                            spath += appb
                            count += 1
                            ipath += 1

                    if verbose > 20:
                        print("-"*60)
                        print(f'Total {spath:.6f}')
                        print()

            print_result(name = f'-<<{op_a};{op_b}>>', value = f'{-vpathT:.6f}')
            lineal_responses[f'<<{op_a};{op_b}>>'] = vpathT

    if verbose > 10:
        name = f"Lineal Response"
        time_object.add_name_delta_time(name = name,  delta_time=(time() - start))

    return lineal_responses