from libr import *

def calculate_quadratic_response(operator_a: list = None, operator_b: list = None,
                            operator_c: list = None,
                            n_mo_occ: int = None, n_mo_virt: int = None,
                            principal_propagator_a: np.array = None,
                            principal_propagator_b: np.array = None,
                            gpvs: dict = None,
                            time_object: drv_time = None,
                            verbose: int = 0):
    """
    Calculate of the path and total value of stactic quadratic response

    <<A; B, C>> = sum_k sum_n [<0|A|k>PP_a^-1(<k|B|n> - d_{kn} <0|B|0>)PP_b^-1<n|C|0>
                                + <0|C|n>PP_a^-1(<n|B|k> - d_{kn} <0|B|0>)PP_b^-1<k|A|0>
                                + <0|C|k>PP_a^-1(<k|A|n> - d_{kn} <0|A|0>)PP_b^-1<n|B|0>]
    Args:
    ----
    operator_a (list): Name of the first operator
    operator_b (list): Name of the second operator
    operator_c (list): Name of the third operator
    n_mo_occ (int): Ocuppied molecular orbitals
    n_mo_virt (int): Virtual molecular orbitals
    principal_propagator_a (np.array): Inverse of the first principal propagator
    principal_propagator_b (np.array): Inverse of the second principal propagator
    time_object (drv_time): Manage time calculation
    verbose (int): Print level
    """
    rotations: int = n_mo_occ*n_mo_virt

    start = time()

    for index_a, op_a in enumerate(operator_a):
        for index_b, op_b in enumerate(operator_b):
            if index_a > index_b and op_a.split()[0] == op_b.split()[0]:
                continue
            for index_c, op_c in enumerate(operator_c):
                if index_b > index_c and op_b.split()[0] == op_c.split()[0]:
                    continue

                if verbose > 20:
                    print_subtitle(name = f"PATHS: <<{op_a};{op_b},{op_c}>>")

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
                                    gpvs[op_a][a+i*n_mo_virt]\
                                    *gpvs[op_b][b+j*n_mo_virt]
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

            print_result(name = f'-<<{op_a};{op_b}>>', value = f'{-vpathT:.6f}')

    if verbose > 10:
        name = f"Lineal Response"
        time_object.add_name_delta_time(name = name,  delta_time=(time() - start))