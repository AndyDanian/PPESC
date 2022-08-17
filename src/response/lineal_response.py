from libr import *
from numba import njit

# Equivalente a indicar
# @jit(nopython=True)
# @njit
def get_lineal_response(n_mo_occ: int = None, n_mo_virt: int = None,
                        n_rotation: int = None,
                        opa: np.array = None, opb: np.array = None,
                        pp: np.array = None, verbose: int = 0):
    ipath: int = 0
    vpathT: float = 0.0E+0

    i: int = 0
    a: int = 0
    for ia in range(n_rotation):
        s = a + n_mo_occ

        count = 0
        spath = 0.0E+0
        j: int = 0
        b: int = 0
        for jb in range(n_rotation):
            t = b + n_mo_occ

            # -0.5*(<i|A|a>P_{ia,jb}^{-1}<b|B|j>
            # +
            # <i|B|b>P_{jb,ia}^{-1}<a|A|j>)
            appb = -0.5*(
                    opa[ia]*
                    pp[ia,jb]*
                    opb[jb]
                +
                    opb[jb+n_rotation]*
                    pp[jb,ia]*
                    opa[ia+n_rotation]
                )

            vpathT += appb
            if verbose > 20 and count == 0:
                print()
                print(f" # ".center(8),f"i".center(8),
                    f"s".center(8),
                    f"t".center(8),f"j".center(8))
            if verbose > 20 and abs(appb) > 0.1: #improve threshild to print value according its amount and basis set size
                print(f"{count + 1}".center(8),f"{i + 1}".center(8),
                    f"{s + 1}".center(8),f"{t + 1}".center(8),
                    f"{j + 1}".center(8),appb
                    )

            spath += appb
            count += 1
            ipath += 1

            b += 1
            if b == n_mo_virt:
                b  = 0
                j += 1

        if verbose > 20:
            print("-"*60)
            print("Total ",spath)
            print()
        a += 1
        if a == n_mo_virt:
            a  = 0
            i += 1

    return vpathT

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

            vpathT = get_lineal_response(n_mo_occ=n_mo_occ,
                                        n_mo_virt=n_mo_virt,
                                        n_rotation = n_rotation,
                                        opa = np.array(gpvs[op_a]),
                                        opb = np.array(gpvs[op_b]),
                                        pp = principal_propagator,
                                        verbose=verbose)

            print_result(name = f'-<<{op_a};{op_b}>>', value = f'{-vpathT:.6f}')
            lineal_responses[f'<<{op_a};{op_b}>>'] = vpathT

    if verbose > 10:
        name = f"Lineal Response"
        time_object.add_name_delta_time(name = name,  delta_time=(time() - start))

    return lineal_responses