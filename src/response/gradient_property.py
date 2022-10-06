from libr import *


def print_gradient_property_vector(
    output: str,
    gpvs: dict[str, list[float]],
) -> None:
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
            f.write(f" GPV  {name}".center(70) + "\n")
            f.write(("-" * 20).center(70) + "\n")
            nrotations: int = int(len(gpv) / 2)

            if max(gpv) > 9999 or min(gpv) < -9999:
                formate: str = "{:.6e}"
            else:
                formate = "{:.6f}"

            for i in range(nrotations):
                f.write(
                    f"{i+1}".center(10)
                    + formate.format(gpv[i]).rjust(20)
                    + formate.format(gpv[i + nrotations]).rjust(20)
                    + "\n"
                )
            f.write("\n")


def print_sot(
    output: str, name: str, arrays: dict[str, list[float]], n_mo: int
) -> None:
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
                f.write(f" <i|{label}|j> ".center(70) + "\n")
            else:
                f.write(f" <a|{label}|b> ".center(70) + "\n")
            f.write(("-" * 20).center(70) + "\n")

            if max(values) > 9999.0 or min(values) < -9999.0:
                formate: str = "{:.6e}"
            else:
                formate = "{:.6f}"

            for i in range(n_mo):
                for j in range(n_mo):
                    f.write(
                        f"{i+1}".center(10)
                        + " "
                        + f"{j+1}".center(10)
                        + " "
                        + formate.format(values[j + i * n_mo]).rjust(16)
                        + "\n"
                    )
            f.write("\n")


def gradient_property_vector_rpa(
    io: scratch,
    time_object: drv_time,
    wf: wave_function,
    list_1b_integrals: list,
    quadratic: bool = False,
    verbose: int = 0,
) -> tuple[dict[str, list[float]], dict[str, list[float]], dict[str, list[float]]]:
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
    start: float = time()

    # molecular integrals
    n_mo_occ: int = wf.mo_occ
    n_mo_virt: int = wf.mo_virt
    mo_coeff_T: np.ndarray = np.array(wf.mo_coefficients)

    gpvs: dict[str, list[float]] = {}
    # Quadratic
    # avs: dict = {}
    mo_occupied: dict[str, list[float]] = {}
    mo_virtuals: dict[str, list[float]] = {}
    #
    for name in list_1b_integrals:
        mo_integral = np.matmul(
            mo_coeff_T,
            np.matmul(
                io.binary(file=io._hermite_ao1b_binary, label=name, io="r"),
                mo_coeff_T.T,
            ),
        )
        if quadratic:
            # avs[name] = 2.0*sum([mo_integral[i][i] for i in range(n_mo_occ)])
            # <i|O|j>
            mo_occupied[name] = [
                mo_integral[a][b] for a in range(n_mo_occ) for b in range(n_mo_occ)
            ]
            # <a|O|b>
            mo_virtuals[name] = [
                mo_integral[a + n_mo_occ][b + n_mo_occ]
                for a in range(n_mo_virt)
                for b in range(n_mo_virt)
            ]
        # gradient porperty vector
        #! Python
        # <i|O|a>
        gpvs[name] = [
            2.0 * mo_integral[a + n_mo_occ][i]
            for i in range(n_mo_occ)
            for a in range(n_mo_virt)
        ]
        # <b|O|j>
        gpvs[name] += [
            -2.0 * mo_integral[i][a + n_mo_occ]
            for i in range(n_mo_occ)
            for a in range(n_mo_virt)
        ]
        if verbose > 20:
            time_object.add_name_delta_time(
                name=f"Build GPV {name}", delta_time=(time() - start)
            )

    if verbose > 30:
        print_gradient_property_vector(output=io._output_path, gpvs=gpvs)
        if quadratic:
            # print_avs(avs = avs)
            print_sot(
                output=io._output_path,
                name="occupied",
                arrays=mo_occupied,
                n_mo=n_mo_occ,
            )
            print_sot(
                output=io._output_path,
                name="virtual",
                arrays=mo_virtuals,
                n_mo=n_mo_virt,
            )

    return mo_occupied, mo_virtuals, gpvs


# *******************************************
def drv_gradient_property_vector(
    io: scratch,
    driver_time: drv_time,
    wf: wave_function,
    list_1b_integrals: list,
    quadratic: bool = False,
    verbose: int = 0,
) -> tuple[dict[str, list[float]], dict[str, list[float]], dict[str, list[float]]]:
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
    # avs: dict = {}
    # temp_avs,
    mo_occupied, mo_virtuals, gpvs = gradient_property_vector_rpa(
        io=io,
        time_object=driver_time,
        wf=wf,
        list_1b_integrals=list_1b_integrals,
        quadratic=quadratic,
        verbose=verbose,
    )
    return mo_occupied, mo_virtuals, gpvs
