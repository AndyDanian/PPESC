from libr import *

def calculate_lineal_reponse(io: scratch = None,
                             time_object: drv_time = None,
                             operator_a: list = None, 
                             operator_b: list = None,
                             pp_multiplicity: np.array = None,
                             gpvs: dict = None,
                             n_mo_occ: int = None,
                             n_mo_virt: int = None,
                             verbose: int = 0):
    """
    Calculate of the path and total value of stactic lineal response

    Args:
    ----
    io (object:scratch): Driver to driver the output and binary files
    time_object (drv_time): Manage time calculation
    operator_a (list): Name of the operator in the left
    operator_b (list): Name of the operator in the right
    pp_multiplicity (strin): Multiplicty of principal propagator
    gpvs (dict): Gradient properties vectors
    n_mo_occ (int): Ocuppied molecular orbitals
    n_mo_virt (int): Virtual molecular orbitals
    verbose (int): Print level
    """
    start = time()

    lineal_responses: dict = {}
    for index_a, op_a in enumerate(operator_a):
        for index_b, op_b in enumerate(operator_b):
            if index_a > index_b and op_a.split()[0] == op_b.split()[0]:
                continue

            vpathT = lineal_sum(# A operator
                                gpva=gpvs[op_a],
                                # B operator
                                gpvb=gpvs[op_b],
                                pp = np.asfortranarray(
                                                        io.binary(file = io._principal_propagator,
                                                                  io = "r",
                                                                  label = pp_multiplicity
                                                                )
                                                        ),
                                # Other
                                verbose=verbose,
                                nocc=n_mo_occ,nvir=n_mo_virt
                                )
            io.write_output(information = f'-<<{op_a};{op_b}>> = {-vpathT:.6f}', type = 1, title_type = 2)
            io.write_output("\n")

            lineal_responses[f'<<{op_a};{op_b}>>'] = vpathT

    if verbose > 10:
        name = f"Lineal Response"
        time_object.add_name_delta_time(name = name,  delta_time=(time() - start))

    return lineal_responses