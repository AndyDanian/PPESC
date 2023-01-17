from libr import *


def calculate_quadratic_response(
    io: scratch,
    time_object: drv_time,
    operator_a: list,
    operator_b: list,
    operator_c: list,
    pp_multiplicity_a: str,
    pp_multiplicity_b: str,
    pp_multiplicity_c: str,
    gpvs: dict,
    mo_occupied: dict,
    mo_virtuals: dict,
    # principal_propagator_a: np.array = None,
    # principal_propagator_b: np.array = None,
    # principal_propagator_c: np.array = None,
    n_mo_occ: int,
    n_mo_virt: int,
    verbose: int = 0,
) -> dict[str, float]:
    """
    Calculate of the path and total value of stactic quadratic response

    <<A; B, C>> = sum_k sum_n [<0|A|k>PP_a^-1(<k|B|n> - d_{kn} <0|B|0>)PP_b^-1<n|C|0>
                                + <0|C|n>PP_a^-1(<n|B|k> - d_{kn} <0|B|0>)PP_b^-1<k|A|0>
                                + <0|C|k>PP_a^-1(<k|A|n> - d_{kn} <0|A|0>)PP_b^-1<n|B|0>
                                + Pertubation]
    Args:
    ----
    io (object:scratch): Driver to driver the output and binary files
    time_object (drv_time): Manage time calculation
    operator_a (list): Name of the first operator
    operator_b (list): Name of the second operator
    operator_c (list): Name of the third operator
    principal_propagator_a (str): Inverse of the first principal propagator save in binary file
    principal_propagator_b (str): Inverse of the second principal propagator save in binary file
    principal_propagator_b (str): Inverse of the third principal propagator save in binary file
    avs (dict): dictionary with average value
    gpvs (dict): Gradient properties vectors
    mo_occupied (dict): Values in molecular orbitals of the occupied orbitals
    mo_virtuals (dict): Values in molecular orbitals of the virtuals orbitals
    n_mo_occ (int): Ocuppied molecular orbitals
    n_mo_virt (int): Virtual molecular orbitals
    verbose (int): Print level
    """
    start: float = time()

    quadratic_responses: dict = {}
    for index_a, op_a in enumerate(operator_a):
        for index_b, op_b in enumerate(operator_b):
            if index_a > index_b and op_a.split()[0] == op_b.split()[0]:
                continue
            for index_c, op_c in enumerate(operator_c):
                if index_b > index_c and op_b.split()[0] == op_c.split()[0]:
                    continue

                vpathT = quadratic_sum(
                    # A operator
                    gpva=gpvs[op_a],
                    iaj=mo_occupied[op_a],
                    sat=mo_virtuals[op_a],
                    ppa=np.asfortranarray(
                        io.binary(
                            file=io._principal_propagator,
                            io="r",
                            label=pp_multiplicity_a,
                        )
                    ),
                    # B operator
                    gpvb=gpvs[op_b],
                    ibj=mo_occupied[op_b],
                    sbt=mo_virtuals[op_b],
                    ppb=np.asfortranarray(
                        io.binary(
                            file=io._principal_propagator,
                            io="r",
                            label=pp_multiplicity_b,
                        )
                    ),
                    # C operator
                    gpvc=gpvs[op_c],
                    icj=mo_occupied[op_c],
                    sct=mo_virtuals[op_c],
                    ppc=np.asfortranarray(
                        io.binary(
                            file=io._principal_propagator,
                            io="r",
                            label=pp_multiplicity_c,
                        )
                    ),
                    # Other
                    verbose=verbose,
                    nocc=n_mo_occ,
                    nvir=n_mo_virt,
                )

                # Following Lehman representation [39]
                # [39] H. Lehman, Nuovo Cim. 11 (1954) 342.
                # Equation (2.43) in Compt. Phys. Reports. 1984, 2, 33
                # <<A;B>>_w = Sum_{n!=0} [<0|A|n><n|B|0> + <0|B|n><n|A|0>](E_0-E_n)^{-1}
                io.write_output(
                    information=f"<<{op_a};{op_b},{op_c}>> = {vpathT:.6e}",
                    type=1,
                    title_type=2,
                )
                io.write_output("\n")
                quadratic_responses[f"<<{op_a};{op_b},{op_c}>>"] = vpathT

    if verbose > 10:
        name = f"Quadratic Response"
        time_object.add_name_delta_time(name=name, delta_time=(time() - start))

    return quadratic_responses
