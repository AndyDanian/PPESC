from libr import *

def print_gradient_property_vector(output: object = None, gpvs: list = None,):
    """
    Print gradient ptoperty vector

    Args:
    ----
        output (object): Output file
        gpvs (list): Values of gradient property vector
        multiplicity (str, int): Multiplicty, open/closed system
    """
    with open(output, "a") as f:
        for name, gpv in gpvs.items():
            f.write(f" GPV  {name}".center(70)+"\n")
            f.write(("-"*20).center(70)+"\n")
            nrotations: int = int(len(gpv)/2)
            if max(gpv) > 9999 or min(gpv) < -9999: 
                formate: str = "{:.6e}"
            else:
                formate: str = "{:.6f}"
            for i in range(nrotations):
                f.write(f"{i+1}".center(10)
                        + formate.format(gpv[i]).rjust(20)
                        + formate.format(gpv[i+nrotations]).rjust(20)
                        + "\n")
            f.write("\n")

def print_avs(avs: list = None):
    """
    Print averages values

    Args:
    ----
        avs (dict): Average values
    """
    for name, value in avs.items():
        print_result(name = f" Average Value {name}", value = value)

def print_sot(output: object = None, 
            name: str = None,
            arrays: list = None,
            n_mo: int = None):
    """
    Print virtuales values in molecular orbital base

    Args:
    ----
        output (object): Output file
        name (str): Occupied or Virtuals values
        arrays (dict): values in molecular orbital
        n_mo (int): Number of occupied or virtuals orbitals
    """
    with open(output, "a") as f:
        for label, values in arrays.items():
            if name == "occupied":
                f.write(f" <i|{label}|j> ".center(70)+"\n")
            else:
                f.write(f" <a|{label}|b> ".center(70)+"\n")
            f.write(("-"*20).center(70)+"\n")
            if max(values) > 9999. or min(values) < -9999.:
                formate: str = "{:.6e}"
            else:
                formate: str = "{:.6f}"

            for i in range(n_mo):
                for j in range(n_mo):
                    f.write(f"{i+1}".center(10) + " " +
                            f"{j+1}".center(10) + " " +
                            formate.format(values[j+i*n_mo]).rjust(16) + "\n")
            f.write("\n")
def gradient_property_vector_rpa(io: scratch = None,
                                time_object: drv_time = None,
                                wf: wave_function = None, 
                                list_1b_integrals: list = None,
                                quadratic: bool = None, 
                                verbose: int = 0,
                                ):
    """
    Calculate of gradient property vectos in rpa approximation

    Args:
    ----
        io (object:scratch): Driver to driver the output and binary files
        time_object (drv_time): Object relationed with the manage of the time
        wf (wave_function): Wave Function object
        gauge (list): Gauge coordiante
        verbose (int): Print level.
    """
    start = time()

    # molecular integrals
    n_mo_occ = wf.mo_occ
    n_mo_virt = wf.mo_virt
    mo_coeff_T = np.array(wf.mo_coefficients)

    gpvs: dict = {}
    # Quadratic
    #avs: dict = {}
    mo_occupied: dict = {}
    mo_virtuals: dict = {}
    #
    for name in list_1b_integrals:
        mo_integral = np.matmul(mo_coeff_T,
                                np.matmul(io.binary(file = io._hermite_ao1b_binary,
                                                    label = name,
                                                    io = "r"),
                                mo_coeff_T.T)
                                )
        if quadratic:
            #avs[name] = 2.0*sum([mo_integral[i][i] for i in range(n_mo_occ)])
            # <i|O|j>
            mo_occupied[name] = [mo_integral[a][b] for a in range(n_mo_occ) for b in range(n_mo_occ)]
            # <a|O|b>
            mo_virtuals[name] = [mo_integral[a + n_mo_occ][b + n_mo_occ]
                                for a in range(n_mo_virt) for b in range(n_mo_virt)
                                ]
        # gradient porperty vector
        #! Python 
        # <i|O|a>
        gpvs[name] = [2.0*mo_integral[a + n_mo_occ][i] for i in range(n_mo_occ) for a in range(n_mo_virt)]
        # <b|O|j>
        gpvs[name] += [-2.0*mo_integral[i][a + n_mo_occ]
                        for i in range(n_mo_occ) for a in range(n_mo_virt)]
        if verbose > 20:
            time_object.add_name_delta_time(name = f"Build GPV {name}", delta_time = (time() - start))

    if verbose > 30:
        print_gradient_property_vector(output = io._output_path, gpvs = gpvs)
        if quadratic:
            #print_avs(avs = avs)
            print_sot(output = io._output_path, name = "occupied", arrays = mo_occupied, n_mo = n_mo_occ)
            print_sot(output = io._output_path, name = "virtual", arrays = mo_virtuals, n_mo = n_mo_virt)

    return mo_occupied, mo_virtuals, gpvs

def read_gradient_property_vector_rpa(gpvs: dict = None, property: str = None,
                                    verbose: int = 0):
    """
    Calculate of gradient property vectos in rpa approximation

    Args:
    ----
        gpvs (dict): Gradient property vectors
        property (str): Property label.
        verbose (int): Print level GPV
        verbose_int (int): Print level for atomic integrals

    """

    read_gpvs = {}
    for name, gpv in gpvs.items():
        if name == property:
            read_gpvs[name] = gpv

    if verbose > 30:
        print_gradient_property_vector(gpv = read_gpvs)

    return read_gpvs

#*******************************************
def drv_gradient_property_vector(io: scratch = None,
                                driver_time: drv_time = None,
                                wf: wave_function = None,
                                list_1b_integrals: list = None,
                                quadratic: bool = False,
                                verbose: int = 0
                                ):
    """
    Driver to build gradient property vector

    Args:
    ----
    io (object:scratch): Driver to driver the output and binary files
    driver_time (drv_time): Object relationed with the manage of the time
    wf (wave_function): Wave function object
    properties (list): Label of the property
    average: (bool): Activate the average value when is calculated quadratic
                    response
    gauge (list): Gauge coordinate
    verbose (int): Print level
    verbose_integrals (int): Print level for hermite module
    """
    #avs: dict = {}
    gpvs: dict = {}
    mo_occupied: dict = {}
    mo_virtuals: dict = {}
    #temp_avs,
    mo_occupied, mo_virtuals, gpvs =\
                            gradient_property_vector_rpa(io = io,
                                                        time_object = driver_time,
                                                        wf = wf,
                                                        list_1b_integrals= list_1b_integrals, 
                                                        quadratic = quadratic,
                                                        verbose = verbose, 
                                                        )
    return mo_occupied, mo_virtuals, gpvs
