from libr import *

def get_principal_propagator_lineal_rpa(n_mo_occ: int = None, n_mo_virt: int = None,
                                        moe: np.array = None, coulomb: np.array = None,
                                        exchange: np.array = None, multiplicity: str = None,
                                        tp_inv: int = 0, verbose: int = 0):
    """
    Build principal propagator for lineal response at rpa level

    Args:
    ----
    n_mo_occ (int): Ocuppied molecular orbitals
    n_mo_virt (int): Virtual molecular orbitals
    moe (np.array, 1d): Molecular orbital energies
    coulomb (np.arra, 4d): Coulomb integrals
    exchange (np.array): Exchange integrals
    multiplicity (str): Multiplicity response
    tp_inv (int): Type of inverse: 0/numpy or 1/series
    verbose (int): Print level
    """

    D1: float = 1.0E+0
    k: int = 0
    irow: int = 0
    icol: int = 0
    rotations: int = n_mo_virt*n_mo_occ

    if tp_inv == 1:
        fock_inv: np.array = np.zeros((rotations,rotations),dtype=float)
    w: np.array = np.zeros((rotations,rotations),dtype=float)

    start = time()
    for i in range(n_mo_occ):
        for a in range(n_mo_virt):
            s = a + n_mo_occ
            if tp_inv == 1:
                fock_inv[k,k] = D1/(moe[s]-moe[i])
            k += 1

            for j in range(n_mo_occ):
                for b in range(n_mo_virt):

                    delta_moe = 0.0
                    if i == j and a == b and tp_inv == 0:
                        delta_moe = moe[s]-moe[i]
                    if (isinstance(multiplicity, str) and multiplicity.lower() == "singlet") or\
                        (isinstance(multiplicity, int) and multiplicity == 1):

                        w[irow,icol] =\
                            -delta_moe - exchange[a,j,b,i] + coulomb[a,b,j,i]

                    elif (isinstance(multiplicity, str) and multiplicity.lower() == "triplet") or\
                        (isinstance(multiplicity, int) and multiplicity == 3):  #<ab|ji> + <aj|bi>

                        w[irow,icol] = -delta_moe + coulomb[a,b,j,i] + exchange[a,j,b,i]

                    else:
                        raise ValueError("***ERROR\n\n\
                            Multiplicity can be 1/singlet or 3/triplet")

                    icol = icol + 1
                    if icol > rotations-1 : irow += 1
                    if icol > rotations-1 : icol  = 0

    build_time = time() - start

    start = time()
    if tp_inv == 0:
        pp = np.linalg.inv(w)
    else:
        raise ValueError("***ERROR\n\n\
            Inverse of a matrix only with numpy is implemeted")
    inverse_time = time() - start

    if verbose > 10:
        print(f"Time take to build principal propagator: {build_time}")
        print(f"Time to inverse calculate: {inverse_time}")

    if verbose > 20:
        if tp_inv == 0:
            name = "Principal Propagator Inverse (Numpy)."
        print_triangle_matrix(integral = pp, name = name,  matriz_sym = "square")
    return pp