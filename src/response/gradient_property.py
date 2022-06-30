from libr import *

def print_gradient_property_vector(gpvs: list = None, multiplicity: str or int = None):
    """
    Print gradient ptoperty vector

    Args:
    ----
        gpvs (list): Values of gradient property vector
        multiplicity (str, int): Multiplicty, open/closed system
    """
    for name, gpv in gpvs.items():
        print()
        print(f" GPV  {name}".center(70))
        print(("-"*20).center(70))
        for count, value in enumerate(gpv):
            if multiplicity in ["Triplet", "triplet", 3]:
                print(f"{count+1}".center(10),f"{value:.6f}".rjust(20),f"{-value:.6f}".rjust(20))
            else:
                print(f"{count+1}".center(10),f"{value:.6f}".rjust(20),f"{value:.6f}".rjust(20))


def gradient_property_vector_rpa(wf: wave_function = None, property: str = None,
                                    multiplicity: str or int = None,
                                    time_object: drv_time = None,
                                    verbose: int = 0, verbose_int: int = 0):
    """
    Calculate of gradient property vectos in rpa approximation

    Args:
    ----
        wf (wave_function): Wave Function object
        property (str): Property label.
        multiplicity (str, int): Multiplicty, open/closed system
        verbose (int): Print level.
    """
    start = time()
    all_responses: bool = False
    # atomic integrals
    calculate_integral = eint(wf)
    integrals_1b, symmetries_1b = calculate_integral.integration_onebody(
    integrals_names = [property], output = verbose_int)

    # molecular integrals
    n_mo_occ = wf.mo_occ
    n_mo_virt = wf.mo_virt
    mo_coeff_T = np.array(wf.mo_coefficients) #from molden come mo coefficient transpose

    gpvs = {}
    for name, atomic_int in integrals_1b.items():
        mo_integral = np.matmul(mo_coeff_T,np.matmul(np.array(atomic_int), mo_coeff_T.T))
        # gradient porperty vector
        gpvs[name] = [2.0*mo_integral[i][a + n_mo_occ] for i in range(n_mo_occ) for a in range(n_mo_virt)]

    if len(gpvs) > 1: all_responses = True # To activate responses among all integrals

    if verbose > 10:
        time_object.add_name_delta_time(name = f"Build GPV {property}", delta_time = (time() - start))

    if verbose > 30:
        print_gradient_property_vector(gpvs = gpvs, multiplicity = multiplicity)

    return all_responses, gpvs

def read_gradient_property_vector_rpa(gpvs: dict = None, property: str = None,
                                    multiplicity: str or int = None, verbose: int = 0):
    """
    Calculate of gradient property vectos in rpa approximation

    Args:
    ----
        gpvs (dict): Gradient property vectors
        property (str): Property label.
        multiplicity (str, int): Multiplicty, open/closed system
        verbose (int): Print level GPV
        verbose_int (int): Print level for atomic integrals

    """

    read_gpvs = {}
    for name, gpv in gpvs.items():
        if name == property:
            read_gpvs[name] = gpv

    if verbose > 30:
        print_gradient_property_vector(gpv = read_gpvs, multiplicity = multiplicity)

    return read_gpvs