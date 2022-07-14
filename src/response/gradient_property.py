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
        nrotations: int = int(len(gpv)/2)
        for i in range(nrotations):
            print(f"{i+1}".center(10),f"{gpv[i]:.6f}".rjust(20),f"{gpv[i+nrotations]:.6f}".rjust(20))

def print_avs(avs: list = None):
    """
    Print averages values

    Args:
    ----
        avs (dict): Average values
    """
    for name, value in avs.items():
        print_result(name = f" Average Value {name}", value = value)

def print_virtuals(virtuals: list = None, n_mo_v: int = None):
    """
    Print virtuales values in molecular orbital base

    Args:
    ----
        virtuals (dict): Virtual value in molecular orbital
    """
    for name, values in virtuals.items():
        print_subtitle(name = f" Virtuales Values {name}")
        for i in range(n_mo_v):
            print(*[f"{values[j+i*n_mo_v]:.6f}" for j in range(n_mo_v)],end="\n")
        print()

def gradient_property_vector_rpa(wf: wave_function = None, property: str = None,
                                    multiplicity: str or int = None,
                                    time_object: drv_time = None,
                                    average: bool = None, gaugeo: list = None,
                                    verbose: int = 0, verbose_int: int = 0):
    """
    Calculate of gradient property vectos in rpa approximation

    Args:
    ----
        wf (wave_function): Wave Function object
        property (str): Property label.
        multiplicity (str, int): Multiplicty, open/closed system
        gaugeo (list): Gauge coordiante
        verbose (int): Print level.
    """
    start = time()
    # atomic integrals
    calculate_integral = eint(wf)
    integrals_1b, symmetries_1b = calculate_integral.integration_onebody(
    integrals_names = [property], verbose = verbose_int, gaugeo = gaugeo)

    # molecular integrals
    n_mo_occ = wf.mo_occ
    n_mo_virt = wf.mo_virt
    mo_coeff_T = np.array(wf.mo_coefficients) #from molden come mo coefficient transpose

    gpvs: dict = {}
    # Quadratic
    avs: dict = {}
    mo_occupied: dict = {}
    mo_virtuals: dict = {}
    #
    for name, atomic_int in integrals_1b.items():
        mo_integral = np.matmul(mo_coeff_T,np.matmul(np.array(atomic_int), mo_coeff_T.T))
        if average:
            avs[name] = 2.0*sum([mo_integral[i][i] for i in range(n_mo_occ)])
            mo_occupied[name] = [mo_integral[a][b] for a in range(n_mo_occ) for b in range(n_mo_occ)]
            mo_virtuals[name] = [mo_integral[a + n_mo_occ][b + n_mo_occ] for a in range(n_mo_virt) for b in range(n_mo_virt)]
        # gradient porperty vector
        gpvs[name] = [2.0*mo_integral[i][a + n_mo_occ] for i in range(n_mo_occ) for a in range(n_mo_virt)]
        gpvs[name] += [-value if multiplicity in ["Triplet", "triplet", 3] else value for value in gpvs[name]]

    if verbose > 10:
        time_object.add_name_delta_time(name = f"Build GPV {property}", delta_time = (time() - start))

    if verbose > 30:
        print_gradient_property_vector(gpvs = gpvs, multiplicity = multiplicity)
        if average:
            print_avs(avs = avs)
            print_virtuals(virtuals = mo_occupied, n_mo_v = n_mo_occ)
            print_virtuals(virtuals = mo_virtuals, n_mo_v = n_mo_virt)

    return avs, mo_occupied, mo_virtuals, gpvs

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

#*******************************************
def drv_gradient_property_vector(wf: wave_function = None, properties: list = None,
                                properties_multiplicity: dict = None,
                                gpv_in: dict = None, driver_time: drv_time = None,
                                average: bool = False, gaugeo: list = None,
                                verbose: int = 0, verbose_integrals: int = -1):
    """
    Driver to build gradient property vector

    Args:
    ----
    wf (wave_function): Wave function object
    properties (list): Label of the property
    properties_multiplicity (dict): Multiplicity for each property
    driver_time (drv_time): Object relationed with the manage of the time
    gpv_in (dict): Gradient properties vectors to read
    average: (bool): Activate the average value when is calculated quadratic
                    response
    gaugeo (list): Gauge coordinate
    verbose (int): Print level
    verbose_integrals (int): Print level for hermite module
    """
    gpvs: dict = {}
    avs: dict = {}
    mo_occupied: dict = {}
    mo_virtuals: dict = {}
    for property in properties:
        if not gpv_in or property not in gpv_in.keys():
            temp_avs, temp_mo_occupied, temp_mo_virtuals, temp_gpvs =\
                                gradient_property_vector_rpa(wf = wf,
                                property = property, time_object = driver_time,
                                verbose = verbose, verbose_int = verbose_integrals,
                                multiplicity=properties_multiplicity[property],
                                average = average, gaugeo = gaugeo
                                )
            gpvs.update(temp_gpvs)
            mo_virtuals.update(temp_mo_virtuals)
            mo_occupied.update(temp_mo_occupied)
            avs.update(temp_avs)
        else:
            gpvs = read_gradient_property_vector_rpa(wf = gpv_in, property = property,
                                verbose = verbose,
                                multiplicity=properties_multiplicity[property])
            raise ValueError("Falta implementar leer los valores average y average virtuales")
    return avs, mo_occupied, mo_virtuals, gpvs