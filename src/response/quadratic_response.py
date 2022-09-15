from libr import *

def calculate_quadratic_response(operator_a: list = None, operator_b: list = None,
                            operator_c: list = None,
                            n_mo_occ: int = None, n_mo_virt: int = None,
                            principal_propagator_a: np.array = None,
                            principal_propagator_b: np.array = None,
                            principal_propagator_c: np.array = None,
                            avs: dict = None, mo_occupied: dict = None,
                            mo_virtuals: dict = None ,gpvs: dict = None,
                            time_object: drv_time = None,
                            verbose: int = 0):
    """
    Calculate of the path and total value of stactic quadratic response

    <<A; B, C>> = sum_k sum_n [<0|A|k>PP_a^-1(<k|B|n> - d_{kn} <0|B|0>)PP_b^-1<n|C|0>
                                + <0|C|n>PP_a^-1(<n|B|k> - d_{kn} <0|B|0>)PP_b^-1<k|A|0>
                                + <0|C|k>PP_a^-1(<k|A|n> - d_{kn} <0|A|0>)PP_b^-1<n|B|0>
                                + Pertubation]
    Args:
    ----
    operator_a (list): Name of the first operator
    operator_b (list): Name of the second operator
    operator_c (list): Name of the third operator
    n_mo_occ (int): Ocuppied molecular orbitals
    n_mo_virt (int): Virtual molecular orbitals
    principal_propagator_a (np.array): Inverse of the first principal propagator
    principal_propagator_b (np.array): Inverse of the second principal propagator
    principal_propagator_b (np.array): Inverse of the third principal propagator
    avs (dict): dictionary with average value
    mo_occupied (dict): Values in molecular orbitals of the occupied orbitals
    mo_virtuals (dict): Values in molecular orbitals of the virtuals orbitals
    time_object (drv_time): Manage time calculation
    verbose (int): Print level
    """
    start: float = time()

    quadratic_responses: dict = {}
    nrot: int = n_mo_occ*n_mo_virt
    for index_a, op_a in enumerate(operator_a):
        for index_b, op_b in enumerate(operator_b):
            if index_a > index_b and op_a.split()[0] == op_b.split()[0]:
                continue
            for index_c, op_c in enumerate(operator_c):
                if index_b > index_c and op_b.split()[0] == op_c.split()[0]:
                    continue

                vpathT = quadratic_sum(gpva=gpvs[op_a],gpvb=gpvs[op_b],gpvc=gpvs[op_c],
                                        iaj=mo_occupied[op_a],
                                        ibj=mo_occupied[op_b],
                                        icj=mo_occupied[op_c],
                                        sat=mo_virtuals[op_a],
                                        sbt=mo_virtuals[op_b],
                                        sct=mo_virtuals[op_c],
                                        ppa = np.asfortranarray(principal_propagator_a),
                                        ppb = np.asfortranarray(principal_propagator_b),
                                        ppc = np.asfortranarray(principal_propagator_c),
                                        verbose=verbose,
                                        nocc=n_mo_occ,nvir=n_mo_virt
                                        )

                print_result(name = f'-<<{op_a};{op_b},{op_c}>>', value = f'{-vpathT:.6f}')
                quadratic_responses[f'<<{op_a};{op_b},{op_c}>>'] = vpathT


    if verbose > 10:
        name = f"Quadratic Response"
        time_object.add_name_delta_time(name = name,  delta_time=(time() - start))

    return quadratic_responses