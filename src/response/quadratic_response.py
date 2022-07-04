from libr import *

def calculate_quadratic_response(operator_a: list = None, operator_b: list = None,
                            operator_c: list = None,
                            n_mo_occ: int = None, n_mo_virt: int = None,
                            principal_propagator_a: np.array = None,
                            principal_propagator_b: np.array = None,
                            avs: dict = None, mo_virtuals: dict = None ,gpvs: dict = None,
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
    avs (dict): dictionary with average value
    mo_virtuals (dict): Values in molecular orbitals of the virtuals orbitals
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
                for i in range(n_mo_occ):
                    #  ---- Virtual index ----
                    for a in range(n_mo_virt):
                        s = a + n_mo_occ
                        count = 0
                        spath = 0.0E+0
                        for b in range(n_mo_virt):
                            t = b + n_mo_occ
                            for c in range(n_mo_virt):
                                u = c + n_mo_occ
                                for d in range(n_mo_virt):
                                    v = d + n_mo_occ
                    #  ---- END Virtual ----
                                    for j in range(n_mo_occ):

                                        if b == c:
                                            vavs_b = avs[op_b]
                                            vavs_a = avs[op_a]
                                        else:
                                            vavs_b = 0.0
                                            vavs_a = 0.0
                                        appb1 =\
                                            (gpvs[op_a][a+i*n_mo_virt]
                                            # PP_{ia,jb}^-1
                                            *principal_propagator_a[a+i*n_mo_virt,b+j*n_mo_virt]
                                            *(mo_virtuals[op_b][c+b*n_mo_virt] - vavs_b)
                                            # PP_{ic,jd}^-1
                                            *principal_propagator_b[c+i*n_mo_virt,d+j*n_mo_virt]
                                            *gpvs[op_c][d+j*n_mo_virt])
                                        appb2 =\
                                            (gpvs[op_c][c+i*n_mo_virt]
                                            # PP_{ia,jb}^-1
                                            *principal_propagator_a[c+i*n_mo_virt,d+j*n_mo_virt]
                                            *(mo_virtuals[op_b][d+a*n_mo_virt] - vavs_b)
                                            # PP_{ic,jd}^-1
                                            *principal_propagator_b[a+i*n_mo_virt,b+j*n_mo_virt]
                                            *gpvs[op_a][b+j*n_mo_virt])
                                        appb3 =\
                                            (gpvs[op_b][a+i*n_mo_virt]
                                            # PP_{ia,jb}^-1
                                            *principal_propagator_a[a+i*n_mo_virt,b+j*n_mo_virt]
                                            *(mo_virtuals[op_a][c+b*n_mo_virt] - vavs_a)
                                            # PP_{ic,jd}^-1
                                            *principal_propagator_b[c+i*n_mo_virt,d+j*n_mo_virt]
                                            *gpvs[op_c][d+j*n_mo_virt])
                                        appb = appb1 + appb2 + appb3

                                        vpathT += appb

                                        if verbose > 20 and count == 0:
                                            print(f" # ".center(6),f"i".center(6),
                                                f"s".center(6),
                                                f"t".center(6),f"u".center(6),
                                                f"v".center(6),f"j".center(6),)
                                        if verbose > 20 and abs(appb) > 0.1:
                                            print(f"{count + 1}".center(6),f"{i + 1}".center(6),
                                                f"{s + 1}".center(6),
                                                f"{t + 1}".center(6),f"{u + 1}".center(6),
                                                f"{v + 1}".center(6),f"{j + 1}".center(6),
                                                f"{appb:.6f}".center(16))

                                        spath += appb
                                        count += 1
                                        ipath += 1

                        if verbose > 20:
                            print("-"*60)
                            print(f'Total {spath:.6f}')
                            print()

                print_result(name = f'-<<{op_a};{op_b},{op_c}>>', value = f'{-vpathT:.6f}')

    if verbose > 10:
        name = f"Quadratic Response"
        time_object.add_name_delta_time(name = name,  delta_time=(time() - start))