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

            vpathT = lineal_sum(gpva=gpvs[op_a],gpvb=gpvs[op_b],
                                pp = np.asfortranarray(principal_propagator),
                                verbose=verbose,
                                nocc=n_mo_occ,nvir=n_mo_virt
                                )
            print_result(name = f'-<<{op_a};{op_b}>>', value = f'{-vpathT:.6f}')
            lineal_responses[f'<<{op_a};{op_b}>>'] = vpathT

    if verbose > 10:
        name = f"Lineal Response"
        time_object.add_name_delta_time(name = name,  delta_time=(time() - start))

    return lineal_responses