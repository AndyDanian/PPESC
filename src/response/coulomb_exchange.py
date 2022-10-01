from libr import *

def get_coulomb_exchange_integrals(driver_time: drv_time = None,
                                   io: scratch = None,
                                   wf: wave_function = None,
                                   verbose: int = 0):
    """
    Get exchange and coulomb integrals from atomic tow-body integrals

    Args:
    ----
        drv_time(object:drv_time): Driver time process
        io (object:scratch): Driver to driver the output and binary files
        wf (wave_function): Wave function object
        verbose (int): Print level for Coulomb and Exchange
    Reference:
    ---------
    Program to transform ao to mo i2c
    C N M Pounder
    Theoret. Chim. Acta (Berl.). 1975, 39, 247--253
    """
    start = time()

    # Coulomb and Exchange
    mo_coeff = np.array(wf.mo_coefficients)
    n_mo_occ = wf.mo_occ
    n_mo_virt = wf.mo_virt

    coulomb = np.zeros((n_mo_virt,n_mo_virt,n_mo_occ,n_mo_occ),dtype=float)
    exchange = np.zeros((n_mo_virt,n_mo_occ,n_mo_virt,n_mo_occ),dtype=float)
    coulomb, exchange = ao2mo(moco = np.asfortranarray(mo_coeff),
                             a2i=np.asfortranarray(
                                io.binary(
                                    file = io._hermite_ao2b_binary,
                                    io = "r",
                                    label = "e2pot"
                                )
                             ),
                             nocc=n_mo_occ, nvir=n_mo_virt)

    ## Write in binary file Exchange and Coulomb integrals
    io.binary(file = io._exchange_coulomb,
              dictionary = {"exchange": exchange},
              io = "a")
    io.binary(file = io._exchange_coulomb,
              dictionary = {"coulomb": coulomb},
              io = "a")
    ## Write in the ouput the EXCHCOUL.H5 file
    io.write_output(information = f"Coulomb and Exchange ({io._exchange_coulomb.name})", 
                    type = 3,
                    size_file = io._exchange_coulomb.stat().st_size)

    driver_time.add_name_delta_time(name = "Exchange and Coulomb", delta_time = (time() - start))

    if verbose > 50:

        io.write_output(information = "Coulomb and Exchange", type = 1, title_type = 1)
        io.write_output(information ="J".center(35) + "K".center(35))
        for a in range(n_mo_virt):
            for b in range(n_mo_virt):
                for i in range(n_mo_occ):
                    for j in range(n_mo_occ):
                        io.write_output(
                                    f"{a + n_mo_occ + 1,b + n_mo_occ + 1,i + 1,j + 1}".center(15) +
                                    f"{coulomb[a,b,i,j]:.8f}".center(20) +
                                    f"{a + n_mo_occ + 1,i + 1,b + n_mo_occ + 1,j + 1}".center(15) +
                                    f"{exchange[a,i,b,j]:.8f}".center(20)
                                    )

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
