from libr import *


def calculate_lineal_reponse(
    io: scratch,
    time_object: drv_time,
    operator_a: list[str],
    operator_b: list[str],
    pp_multiplicity: np.ndarray,
    gpvs: dict[str, list[float]],
    n_mo_occ: int,
    n_mo_virt: int,
    verbose: int = 0,
) -> dict[str, float]:
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
    start: float = time()

    lineal_responses: dict[str, float] = {}
    sign: float = 1.0
    for index_a, op_a in enumerate(operator_a):
        for index_b, op_b in enumerate(operator_b):
            if index_a > index_b and op_a.split()[0] == op_b.split()[0]:
                continue

            if pp_multiplicity == 'triplet':
                sign = -1.0

            vpathT = lineal_sum(
                # A operator
                gpva=gpvs[op_a],
                # B operator
                gpvb=gpvs[op_b],
                pp=np.asfortranarray(
                    io.binary(
                        file=io._principal_propagator, io="r", label=pp_multiplicity
                    )
                ),
                # Other
                verbose=verbose,
                nocc=n_mo_occ,
                nvir=n_mo_virt,
            )
            # Lehman representation [39], Equation (2.43) in Compt. Phys. Reports. 1984, 2, 33
            # [39] H. Lehman, Nuovo Cim. 11 (1954) 342.
            # <<A;B>>_w = Sum_{n!=0} [<0|A|n><n|B|0> + <0|B|n><n|A|0>](E_0-E_n)^{-1}
            # WARNING WARNING WARNING WARNING WARNING WARNING
            #! WARNING: IN THE PATH IS AVOID THE SIGN IN TRIPLET RESPOMSE
            io.write_output(
                information=f"<<{op_a};{op_b}>> = {sign*vpathT:.6e}", type=1, title_type=2
            )
            io.write_output("\n")

            lineal_responses[f"<<{op_a};{op_b}>>"] = sign*vpathT

    if verbose > 10:
        name = f"Lineal Response"
        time_object.add_name_delta_time(name=name, delta_time=(time() - start))

    return lineal_responses
