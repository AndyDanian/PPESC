from libr import *

def calculate_lineal_reponse(n_mo_occ: int = None, n_mo_virt: int = None,
                            principal_propagator: np.array = None, gpvs: dict = None,
                            verbose: int = 0):
    """
    Calculate of the path and total value of lienal response

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
    names_gpv: list = [name for name in gpvs.keys()]
    values_gpv: list = [v for v in gpvs.values()]
    n_gpv: int = len(names_gpv)
    for index_gpv_left in range(n_gpv):
        for index_gpv_right in range(index_gpv_left, n_gpv):

            iset = 0
            ipath = 0
            vpathT = 0.0E+0

            irow = 0
            icol = 0
            for i in range(n_mo_occ):
                for a in range(n_mo_virt):
                    s = a + n_mo_occ
            #!
                    iset += 1
                    count = 0
                    spath = 0.0E+0
            #!
                    for j in range(n_mo_occ):
                        for b in range(n_mo_virt):
                            t = b + n_mo_occ
            #!<i|FC(C)|s><t|FC(H)|j>
            #!2 D2 por ias o tbj, es necesario para reproducir valores del DALTON
                            appb =\
                                values_gpv[index_gpv_left][a+i*n_mo_virt]\
                                    *values_gpv[index_gpv_right][b+j*n_mo_virt]
            #!<i|FC(C)|s><t|FC(H)|j> ESTNEi
                            appb = appb*principal_propagator[irow,icol]
            #!Sum path
                            vpath = 0.0E+0
                            vpath = vpath + appb
                            vpathT = vpathT + vpath

                            if verbose > 4 and count == 0: print('   #    t    i    s    j  ')
                            if verbose > 4 and abs(vpath) > 0.1: \
                                    print(count+1,t+1,i+1,s+1,j+1,vpath)
            #!
                            spath = spath + vpath
                            count += 1
            #!
                            if verbose > 2 and ipath/(n_mo_virt*n_mo_virt) == iset :
                                print('Paths Total',spath)
                                print()
            #!
                            ipath += 1

            print()
            print('************************************')
            # label1 = '-<<' + label[itp].replace('\n','') + ';' + label[itp+1].replace('\n','') + '>>'
            # print(label1)
            print(f'RPA         {vpathT:.6f}')
            # print(f'PZOA        {tpzoa:.6f}')
            # dif = tpzoa-vpathT
            # print(f'Cont. I2C   {dif:.6f}')
            print('***********************************')
