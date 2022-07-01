from libr import *

def get_coulomb_exchange_integrals(wf: wave_function = None,
                                    time_object: drv_time = None,
                                    verbose: int = 0,
                                    verbose_int: int = 0):
    """
    Get exchange and coulomb integrals from atomic tow-body integrals

    Args:
    ----
        wf (wave_function): Wave function object
        verbose (int): Print level for Coulomb and Exchange
        verbose_int (int): Print level for atomic integrals
    Reference:
    ---------
    Program to transform ao to mo i2c
    C N M Pounder
    Theoret. Chim. Acta (Berl.). 1975, 39, 247--253
    """
    start = time()
    # atomic integrals
    calculate_integral = eint(wf)
    at2in: np.array = np.array(calculate_integral.integration_twobody(
        integrals_names = ["e2pot"], verbose = verbose_int,
    )["e2pot"])

    # Coulomb and Exchange
    mo_coeff = np.array(wf.mo_coefficients)
    nprim = wf.primitives_number
    n_mo_occ = wf.mo_occ
    n_mo_virt = wf.mo_virt

    coulomb = np.zeros((n_mo_virt,n_mo_virt,n_mo_occ,n_mo_occ),dtype=float)
    exchange = np.zeros((n_mo_virt,n_mo_occ,n_mo_virt,n_mo_occ),dtype=float)
    for a in range(n_mo_virt):
        s = a + n_mo_occ

        #Transform first index
        f = np.zeros((nprim,nprim,nprim),dtype=float)
        for i in range(nprim):
            for j in range(nprim):
                for k in range(nprim):
                    for l in range(nprim):
                        f[i,j,k] += mo_coeff[s,l]*at2in[i,j,k,l]

        #Coulomb
        for b in range(a, n_mo_virt):
            t = b + n_mo_occ

            # Transform second index
            g = np.zeros((nprim,nprim),dtype=float)
            for i in range(nprim):
                for l in range(nprim):
                    for k in range(nprim):
                        g[i,l] += mo_coeff[t,k]*f[i,l,k]

            for j in range(n_mo_occ):

                #Transform third index
                h = np.zeros((nprim),dtype=float)
                for i in range(nprim):
                    for k in range(nprim):
                            h[i] += mo_coeff[j,k]*g[k,i]

                #Transform four-th index to get Coulomb integrals
                for i in range(j, n_mo_occ):
                    for k in range(nprim):
                        coulomb[a,b,j,i] += mo_coeff[i,k]*h[k]
                    if b >= a:
                        coulomb[b,a,i,j] = coulomb[a,b,j,i]
                        coulomb[a,b,i,j] = coulomb[a,b,j,i]
                        coulomb[b,a,j,i] = coulomb[a,b,j,i]

        #Exchange
        for j in range(n_mo_occ):
            # Transform second index
            g = np.zeros((nprim,nprim),dtype=float)
            for i in range(nprim):
                for l in range(nprim):
                    for k in range(nprim):
                        g[i,l] += mo_coeff[j,k]*f[i,l,k]

            for b in range(a, n_mo_virt):
                t  = b + n_mo_occ
                # Transform third index
                h = np.zeros((nprim),dtype=float)
                for i in range(nprim):
                    for k in range(nprim):
                        h[i] += mo_coeff[t,k]*g[k,i]

                # Transform four-th index
                for i in range(n_mo_occ):
                    for k in range(nprim):
                        exchange[a,j,b,i] += mo_coeff[i,k]*h[k]
                    if b > a:
                        exchange[b,i,a,j] = exchange[a,j,b,i]

    if verbose > 10:
        time_object.add_name_delta_time(name = f"Coulomb and Exchange", delta_time = (time() - start))

    if verbose > 50:

        print("Coulomb and Exchange".center(70))
        print(("-"*30).center(70))
        print()
        print("J".center(35),"K".center(35))
        for a in range(n_mo_virt):
            for b in range(n_mo_virt):
                for i in range(n_mo_occ):
                    for j in range(n_mo_occ):
                        print(f"{a + n_mo_occ + 1,b + n_mo_occ + 1,i + 1,j + 1}".center(15),
                            f"{coulomb[a,b,i,j]:.8f}".center(20),
                            f"{a + n_mo_occ + 1,i + 1,b + n_mo_occ + 1,j + 1}".center(15),
                            f"{exchange[a,i,b,j]:.8f}".center(20)
                        )

    return coulomb, exchange


def ao_2_mo(wf: wave_function = None, verbose: int = 0):
    """
    transform all atomic two body integrals to molecular orbitals

    Args:
    ----
        wf (wave_function): Wave function object

    Reference:
    ---------
    Program to transform ao to mo i2c
    C N M Pounder
    Theoret. Chim. Acta (Berl.). 1975, 39, 247--253
    """

    # atomic integrals
    calculate_integral = eint(wf)
    at2in: np.array = np.array(calculate_integral.integration_twobody(
        integrals_names = ["e2pot"], verbose = verbose,
    )["e2pot"])

    # Coulomb and Exchange
    mo_coeff = np.array(wf.mo_coefficients)
    nprim = wf.primitives_number

    mo2i = np.zeros((nprim,nprim,nprim,nprim),dtype=float)
    count = 0
    for a in range(nprim):

        #Transform first index
        f = np.zeros((nprim,nprim,nprim),dtype=float)
        g = np.zeros((nprim,nprim,nprim),dtype=float)
        for i in range(nprim):
            for j in range(nprim):
                for k in range(nprim):
                    for l in range(nprim):
                        f[i,j,k] += mo_coeff[a,l]*at2in[i,j,k,l]

        # Transform second index
        for i in range(nprim):
            for j in range(nprim):
                for b in range(a, nprim):
                    for k in range(nprim):
                        g[b,i,j] += mo_coeff[b,k]*f[i,j,k]

        #Transform third index
        f = np.zeros((nprim,nprim,nprim),dtype=float)
        for i in range(nprim):
            for b in range(a, nprim):
                for j in range(a, nprim):
                    for k in range(nprim):
                        f[b,j,i] += mo_coeff[j,k]*g[b,k,i]

        #Transform four-th index
        for j in range(a, nprim):
            for b in range(a, nprim):
                for i in range(a, nprim):
                    if abs(mo2i[a,b,j,i]) > 1.0E-16:
                        continue
                    for k in range(nprim):
                        mo2i[a,b,j,i] += mo_coeff[i,k]*f[b,j,k]
                    count += 1
                    if b > a and i == j:
                        mo2i[b,a,j,i] = mo2i[a,b,j,i]
                        mo2i[j,i,a,b] = mo2i[a,b,j,i]
                        mo2i[j,i,b,a] = mo2i[a,b,j,i]
                    elif i > j and a == b:
                        mo2i[a,b,i,j] = mo2i[a,b,j,i]
                        mo2i[j,i,a,b] = mo2i[a,b,j,i]
                        mo2i[i,j,a,b] = mo2i[a,b,j,i]
                    elif b > a and i > j:
                        mo2i[b,a,j,i] = mo2i[a,b,j,i]
                        mo2i[a,b,i,j] = mo2i[a,b,j,i]
                        mo2i[b,a,i,j] = mo2i[a,b,j,i]
                        if b != i or b != j or a != i or a != j:
                            mo2i[j,i,b,a] = mo2i[a,b,j,i]
                            mo2i[i,j,a,b] = mo2i[a,b,j,i]
                            mo2i[i,j,b,a] = mo2i[a,b,j,i]

    if verbose > 50:

        print("Molecular Obrital Two Body Integrals".center(70))
        print(("-"*30).center(70))
        print(" Integrals Number Transfromed : ",count)
        for i in range(nprim):
            for j in range(nprim):
                for k in range(nprim):
                    for l in range(nprim):
                        print(f"{i + 1,j + 1,k + 1,l + 1}".center(15),
                            f"{mo2i[i,j,k,l]:.8f}".center(20)
                        )
