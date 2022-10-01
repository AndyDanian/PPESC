from libl import *

def correction_to_calculate(lresc_amounts: list = None, tensor: bool = False):
    """
    Turn on or off lresc_amounts corrections
    Args:
        paramagneitc (list): Paramagnetic amount to calculate
    """
    triplet_lineal_amount: bool = False
    singlet_lineal_amount: bool = False
    triplet_quadratic_amount: bool = False
    singlet_quadratic_amount: bool = False
    para_nr: bool = False
    dia_nr: bool = False
    dia_lineal: bool = False
    diaavs: bool = False

    sddxy: bool = False
    sddxz: bool = False
    sddyx: bool = False
    sddyz: bool = False
    sddzx: bool = False
    sddzy: bool = False
    lsdsoxy: bool = False
    lsdsoxz: bool = False
    lsdsoyx: bool = False
    lsdsoyz: bool = False
    lsdsozx: bool = False
    lsdsozy: bool = False
    sdbsoxy: bool = False
    sdbsoxz: bool = False
    sdbsoyx: bool = False
    sdbsoyz: bool = False
    sdbsozx: bool = False
    sdbsozy: bool = False
    if tensor:
        sddxy: bool = True
        sddxz: bool = True
        sddyx: bool = True
        sddyz: bool = True
        sddzx: bool = True
        sddzy: bool = True
        lsdsoxy: bool = True
        lsdsoxz: bool = True
        lsdsoyx: bool = True
        lsdsoyz: bool = True
        lsdsozx: bool = True
        lsdsozy: bool = True
        sdbsoxy: bool = True
        sdbsoxz: bool = True
        sdbsoyx: bool = True
        sdbsoyz: bool = True
        sdbsozx: bool = True
        sdbsozy: bool = True


    if not lresc_amounts:
        return {
        "paranr": True,
        # Lineal -- Triplet
        "fclap": True, "sddxx": True, "sddyy": True, "sddzz": True,
        "sddxy": sddxy, "sddxz": sddxz, "sddyx": sddyx,
        "sddyz": sddyz, "sddzx": sddzx, "sddzy": sddzy,
        "fcbso": True, "sdbsoxx": True, "sdbsoyy": True,
        "sdbsozz": True, "sdbsoxy": sdbsoxy, "sdbsoxz": sdbsoxz,
        "sdbsoyx": sdbsoyx, "sdbsoyz": sdbsoyz, "sdbsozx": sdbsozx,
        "sdbsozy": sdbsozy,
        # Lineal -- Singlet
        "lpsokin": True, "lkinpso": True,
        # Quadratic -- Triplet
        "lfcso" : True,
        "lsdsoxx": True, "lsdsoyy": True, "lsdsozz": True,
        "lsdsoxy": lsdsoxy, "lsdsoxz": lsdsoxz, "lsdsoyx": lsdsoyx,
        "lsdsoyz": lsdsoyz, "lsdsozx": lsdsozx, "lsdsozy": lsdsozy,
        # Quadratic -- Singlet
        "lpsodw": True, "lpsomv": True,
        # Diamagnetic
        # NR
        "dianr": True,
        # Lineal singlet responses
        "a2mv": True, "a2dw": True,
        # Averages
        "fc": True, "psooz": True, "dnske": True
        }, {"para_nr": True, "triplet_lineal_amount": True, "singlet_lineal_amount": True,
        "triplet_quadratic_amount": True, "singlet_quadratic_amount": True,
        "dia_singlet_lineal": True}, {"dia_nr": True, "dia_avs": True}
    else:
        if "paranr" in lresc_amounts or "nr" in lresc_amounts: para_nr = True
        if "dianr" in lresc_amounts or "nr" in lresc_amounts: dia_nr = True
        if "diaavs" in lresc_amounts or "avs" in lresc_amounts: diaavs = True
        if ("fclap" in lresc_amounts or "sdlap" in lresc_amounts or
            "fcbso" in lresc_amounts or "sdbso" in lresc_amounts): triplet_lineal_amount = True
        if "lpsokin" in lresc_amounts or "lkinpso" in lresc_amounts: singlet_lineal_amount = True
        if ("lfcso" in lresc_amounts or "lsdsoxx" in lresc_amounts or
            "lsdsoyy" in lresc_amounts or "lsdsozz" in lresc_amounts): triplet_quadratic_amount = True
        if ("lpsodw" in lresc_amounts or "lpsomv" in lresc_amounts): singlet_quadratic_amount = True
        if ("a2mv" in lresc_amounts or "a2dw" in lresc_amounts): dia_lineal = True
        return {
        "paranr": "paranr" in lresc_amounts,
        # Lineal -- Triplet
        "fclap": "fclap" in lresc_amounts, "sddxx": "sddxx" in lresc_amounts,
        "sddyy": "sddyy" in lresc_amounts,"sddzz": "sddzz" in lresc_amounts,
        "sddxy": "sddxy" in lresc_amounts, "sddxz": "sddxz" in lresc_amounts,
        "sddyx": "sddyx" in lresc_amounts, "sddyz": "sddyz" in lresc_amounts,
        "sddzx": "sddzx" in lresc_amounts, "sddzy": "sddzy" in lresc_amounts,
        "fcbso": "fcbso" in lresc_amounts, "sdbsoxx": "sdbsoxx" in lresc_amounts,
        "sdbsoyy": "sdbsoyy" in lresc_amounts, "sdbsozz": "sdbsozz" in lresc_amounts,
        "sdbsoxy": "sdbsoxy" in lresc_amounts, "sdbsoxz": "sdbsoxz" in lresc_amounts,
        "sdbsoyx": "sdbsoyx" in lresc_amounts, "sdbsoyz": "sdbsoyz" in lresc_amounts,
        "sdbsozy": "sdbsozy" in lresc_amounts, "sdbsozx": "sdbsozx" in lresc_amounts,
        # Lineal -- Singlet
        "lpsokin": "lpsokin" in lresc_amounts, "lkinpso": "lkinpso" in lresc_amounts,
        # Quadratic -- Triplet
        "lfcso" : "lfcso" in lresc_amounts,
        "lsdsoxx": "lsdsoxx" in lresc_amounts, "lsdsoyy": "lsdsoyy" in lresc_amounts,
        "lsdsozz": "lsdsozz" in lresc_amounts, "lsdsoxy": "lsdsoxy" in lresc_amounts,
        "lsdsoxz": "lsdsoxz" in lresc_amounts, "lsdsoyx": "lsdsoyx" in lresc_amounts,
        "lsdsoyz": "lsdsoyz" in lresc_amounts, "lsdsozx": "lsdsozx" in lresc_amounts,
        "lsdsozy": "lsdsozy" in lresc_amounts,
        # Quadratic -- Singlet
        "lpsodw": "lpsodw" in lresc_amounts, "lpsomv": "lpsomv" in lresc_amounts,
        # Diamagnetic
        # NR
        "dianr": "dianr" in lresc_amounts,
        # Lineal singlet responses
        "a2mv": "a2mv" in lresc_amounts, "a2dw": "a2dw" in lresc_amounts,
        # Averages
        "fc": "fc" in lresc_amounts, "psooz": "psooz" in lresc_amounts, "dnske": "dnske" in lresc_amounts
        }, {"para_nr": para_nr, "triplet_lineal_amount": triplet_lineal_amount, "singlet_lineal_amount": singlet_lineal_amount,
        "triplet_quadratic_amount": triplet_quadratic_amount, "singlet_quadratic_amount": singlet_quadratic_amount,
        "dia_singlet_lineal": dia_lineal}, {"dia_nr": dia_nr, "dia_avs": diaavs}

def run_shielding_lresc(io: scratch = None,
                        wf: wave_function = None, lresc_amounts: list = None,
                        lresc_consts: dict = None, atom: list = None,
                        principal_propagator_approximation: str = None,
                        driver_time: drv_time = None, verbose: int = 0,
                        tensor: bool = False,
                        verbose_response: int = -1,
                        verbose_average: int = -1,):
    """
    Calculate all lresc_amounts values

    The different amount involve in the LRESC methodology for Shielding
    calculation are split in correction set (dictionaries)

    label_amounts/type_correction=
        {
        para_nr: NR paramagnetic amounts
                {"paranr": paranr}
                            paranr = [
                                    lambda a: ["angmom x", "pso " + str(1 + 3*a)],
                                    lambda a: ["angmom y", "pso " + str(1 + 3*a)],
                                    lambda a: ["angmom z", "pso " + str(1 + 3*a)],
                                    lambda a: ["angmom x", "pso " + str(2 + 3*a)],
                                    lambda a: ["angmom y", "pso " + str(2 + 3*a)],
                                    lambda a: ["angmom z", "pso " + str(2 + 3*a)],
                                    lambda a: ["angmom x", "pso " + str(3 + 3*a)],
                                    lambda a: ["angmom y", "pso " + str(3 + 3*a)],
                                    lambda a: ["angmom z", "pso " + str(3 + 3*a)],
                                    ]
        dia_nr : NR diamagnetic amounts
                {"dianr": dianr}
        triplet_lineal_amount: Correction type triplet response associated with
                                paramagnetic component
                                {"fclap": fclap,
                                "sddxx": sddxx, "sddyy": sddyy, "sddzz": sddzz,
                                ...}
        singlet_lineal_amount: Correction type singlet response associated with
                                paramagnetic component
                                {"lpsokin": lpsokin, "lkinpso": lkinpso}
        triplet_quadratic_amount: Correction type triplet response associated with
                                paramagnetic component and quadratic response
        singlet_quadratic_amount: Correction type singlet response associated with
                                paramagnetic component and quadratic response
        dia_av: Correction type average associated with diamagnetic component
                {"fc": fc, "sd": sd, "psooz": psooz, "dnske": dnske, "pnstcgop": pnstcgop}
        dia_singlet_lineal: Correction type singlet response associated with
                            diamagnetic component
        }

    Also, there are another dictionary, which indicate if the amount will be calculated, it is
    called correction

    Args:
        io (object:scratch): Driver to driver the output and binary files
        wf (wave_function): Wave function object
        lresc_amounts (list): Amounts will calculate
        lresc_consts (dict): Dictionary with constants for differente amounts
        atom (list): Atom index to calculate the shielding
        principal_propagator_approximation (str): Approximation
        tensor (bool): Activate calculation of all tensor components
        driver_time (drv_time): Object to manage the time calculation
        scalar_correction (bool): Activate Dw and Mv correction for energy
        verbose (int): Print level
        verbose_integrals (int): Print level for integral module
    """
    start: float = time()

    corrections, responses_amounts, averages = correction_to_calculate(lresc_amounts,
                                                                    tensor)
    ## The values of next dictionary are in include lresc_parameters.py
    label_amounts: dict = {"para_nr": paramagnetic_nr,
                            "triplet_lineal_amount": triplet_lineal_responses,
                            "singlet_lineal_amount": singlet_lineal_responses,
                            "triplet_quadratic_amount": triplet_quadratic_responses,
                            "singlet_quadratic_amount": singlet_quadratic_responses,
                            "dia_singlet_lineal": dia_singlet_lineal_responses,
                            "dia_nr": diamagnetic_nr,
                            "dia_avs": dia_averages
                            }

    all_averages: dict = {}
    all_responses: dict = {}
    av: average = average(wf)
    moe: list = wf.mo_energies
    lineal_response: response = response(wf = wf, moe = moe)

    delta_average: float = 0.0
    delta_response: float = 0.0
    for a in atom:
        gauge: list = wf.coordinates[a]

        # Avarage calculation
        start_average: float = time()
        if verbose_average < 10: io.activate_write_output = False
        atom_av: dict = {}
        for property, activate in averages.items():
            if activate:
                for name, operators in label_amounts[property].items():
                    if corrections[name]:
                        # Number of tensor components different
                        if name != "fc":
                            tensor_components: int = 9
                        else:
                            tensor_components: int = 1
                        # Diagonal (default) or all components of tensor
                        if not tensor:
                            tensor_step: int = 4
                        else:
                            tensor_step: int = 1
                        temp_av_a: list = [(lresc_consts[name]*
                                                    list(av.calculate_average(
                                                                    property = operators[component](a),
                                                                        gauge = gauge,
                                                                        verbose = verbose_average
                                                                            ).values())[0])
                                        for component in range(0,tensor_components,tensor_step)]
                        if name == "fc":
                            if not tensor:
                                temp_av_a += 2*temp_av_a
                            else:
                                temp: float = temp_av_a[0]
                                temp_av_a: list = []
                                temp_av_a = [temp, 0.0, 0.0, 0.0, temp, 0.0, 0.0, 0.0, temp]
                        atom_av[name] = temp_av_a
        all_averages[a] = atom_av
        if verbose_average < 10: io.activate_write_output = True

        delta_average += time() - start_average
        # End Average calculation

        # Response calculation
        start_response: float = time()
        if verbose_response < 10: io.activate_write_output = False
        atom_responses: dict = {}
        sdlap_components: list =   ["sddxx", "sddxy", "sddxz",
                                    "sddyx", "sddyy", "sddyz",
                                    "sddzx", "sddzy", "sddzz"]
        for responses, activate in responses_amounts.items():
            if activate:
                for label, response_calculation in label_amounts[responses].items():
                    if corrections[label]:
                        # Diagonal (default) or all components of tensor
                        if not tensor:
                            tensor_step: int = 4
                        else:
                            tensor_step: int = 1
                        # Number of tensor components different
                        if (label not in sdlap_components
                            and label != "fclap"):
                            tensor_components: int = 9
                        else:
                            tensor_components: int = 3
                            tensor_step: int = 1
                        temp_responses_a: list = [-lresc_consts[label]
                                                    *list(lineal_response.drv_reponse_calculation(
                                                                                                gauge = gauge,
                                                                principal_propagator_approximation=principal_propagator_approximation,
                                                                            properties = [response_calculation[component](a)],
                                                                                        verbose=verbose_response
                                                                                                ).values())[0]
                                                    for component in range(0,tensor_components,tensor_step)]
                        atom_responses[label] = temp_responses_a
        if not tensor:
            if "sddxx" in atom_responses.keys() and "sddyy" in atom_responses.keys() and "sddzz" in atom_responses.keys():
                atom_responses["sdlap"] =   [
                                            sum(atom_responses["sddxx"]),
                                            sum(atom_responses["sddyy"]),
                                            sum(atom_responses["sddzz"])
                                        ]
                #atom_responses.pop("sddyy")
            if ("sdbsoxx" in atom_responses.keys() and
                "sdbsoyy" in atom_responses.keys() and
                "sdbsozz" in atom_responses.keys()):
                atom_responses["sdbso"] =   [
                                            sum(atom_responses["sdbsoxx"]),
                                            sum(atom_responses["sdbsoyy"]),
                                            sum(atom_responses["sdbsozz"])
                                        ]
            if ("lsdsoxx" in atom_responses.keys() and
                "lsdsoyy" in atom_responses.keys() and
                "lsdsozz" in atom_responses.keys()):
                atom_responses["lsdso"] =   [
                                            sum(atom_responses["lsdsoxx"]),
                                            sum(atom_responses["lsdsoyy"]),
                                            sum(atom_responses["lsdsozz"])
                                        ]
        else:
            if ("sddxx" in atom_responses.keys() and "sddxy" in atom_responses.keys() and "sddxz" in atom_responses.keys()
                and "sddyx" in atom_responses.keys() and "sddyy" in atom_responses.keys() and "sddyz" in atom_responses.keys()
                and "sddzx" in atom_responses.keys() and "sddzy" in atom_responses.keys() and "sddzz" in atom_responses.keys()):
                atom_responses["sdlap"] =   [
                                            sum(atom_responses["sddxx"]),
                                            sum(atom_responses["sddxy"]),
                                            sum(atom_responses["sddxz"]),
                                            sum(atom_responses["sddyx"]),
                                            sum(atom_responses["sddyy"]),
                                            sum(atom_responses["sddyz"]),
                                            sum(atom_responses["sddzx"]),
                                            sum(atom_responses["sddzy"]),
                                            sum(atom_responses["sddzz"]),
                                        ]
            if ("sdbsoxx" in atom_responses.keys() and "sdbsoxy" in atom_responses.keys() and "sdbsoxz" in atom_responses.keys()
                and "sdbsoyx" in atom_responses.keys() and "sdbsoyy" in atom_responses.keys() and "sdbsoyz" in atom_responses.keys()
                and "sdbsozx" in atom_responses.keys() and "sdbsozy" in atom_responses.keys() and "sdbsozz" in atom_responses.keys()):
                atom_responses["sdbso"] =   [
                                            sum(atom_responses["sdbsoxx"]),
                                            sum(atom_responses["sdbsoxy"]),
                                            sum(atom_responses["sdbsoxz"]),
                                            sum(atom_responses["sdbsoyx"]),
                                            sum(atom_responses["sdbsoyy"]),
                                            sum(atom_responses["sdbsoyz"]),
                                            sum(atom_responses["sdbsozx"]),
                                            sum(atom_responses["sdbsozy"]),
                                            sum(atom_responses["sdbsozz"]),
                                        ]
            if ("lsdsoxx" in atom_responses.keys() and "lsdsoxy" in atom_responses.keys() and "lsdsoxz" in atom_responses.keys()
                and "lsdsoyx" in atom_responses.keys() and "lsdsoyy" in atom_responses.keys() and "lsdsoyz" in atom_responses.keys()
                and "lsdsozx" in atom_responses.keys() and "lsdsozy" in atom_responses.keys() and "lsdsozz" in atom_responses.keys()):
                atom_responses["lsdso"] =   [
                                            sum(atom_responses["lsdsoxx"]),
                                            sum(atom_responses["lsdsoxy"]),
                                            sum(atom_responses["lsdsoxz"]),
                                            sum(atom_responses["lsdsoyx"]),
                                            sum(atom_responses["lsdsoyy"]),
                                            sum(atom_responses["lsdsoyz"]),
                                            sum(atom_responses["lsdsozx"]),
                                            sum(atom_responses["lsdsozy"]),
                                            sum(atom_responses["lsdsozz"]),
                                        ]
            if "fclap" in atom_responses.keys():
                temp: float = atom_responses["fclap"]
                atom_responses["fclap"] = [temp[0],0.0,0.0,0.0,temp[1],0.0,0.0,0.0,temp[2]]
        # Remove components of sdlap, sdbso, and lsdso
        all_responses[a] =  {
                                name: atom_responses[name]
                                for name in name_order_responses
                            }
        delta_response += time() - start_response
        if verbose_response < 10: io.activate_write_output = True
        # End Response calculation

    io.activate_write_output = True
    if verbose > 10:
        driver_time.add_name_delta_time(name = "Averages Amounts Calculations", delta_time = delta_average)
        driver_time.add_name_delta_time(name = "Responses Amounts Calculations", delta_time = delta_response)
        driver_time.add_name_delta_time(name = "LRESC Amounts Calculations", delta_time = (time() - start))

    return all_responses, all_averages

def get_shielding_iso_ani(all_responses: dict = None, all_averages: dict = None,
                        tensor: bool = False, ani_axe: str or int = "z"):
    """
    Calculate the isotropic and anisotropic value

    Args:
        all_responses (dict): Label and value of the responses
        all_averages (dict): Label and value of the averages
        tensor (bool): Responses and averages are tensors representation
        ani_axes (str or int): Axes to calculate the anisotropic value

    Returns:
        iso_averages: Isotropic value of averages
        iso_responses: Isotropic value of responses
        ani_averages: Anisotropic value of averages
        ani_responses: Anisotropic value of responses
    """

    sig_x: int = 1.0
    sig_y: int = 1.0
    sig_z: int = -1.0
    if ani_axe == 1 or ani_axe == 0 or ani_axe == "x":
        sig_x: int = -1.0
        sig_y: int = 1.0
        sig_z: int = 1.0
    elif ani_axe == 2 or ani_axe == "y":
        sig_x: int = 1.0
        sig_y: int = -1.0
        sig_z: int = 1.0

    a: int = 0
    b: int = 1
    c: int = 2
    if tensor:
        b: int = 4
        c: int = 8

    isotropic_responses: dict = {}
    anisotropic_responses: dict = {}
    isotropic_averages: dict = {}
    anisotropic_averages: dict = {}
    for atom, amounts in  all_responses.items():
        isotropic: dict = {}
        anisotropic: dict = {}
        for label, responses in amounts.items():
            isotropic[label] = (responses[a] + responses[b] + responses[c])/3.0
            anisotropic[label] = sig_x*responses[a] + sig_y*responses[b] + sig_z*responses[c]
        isotropic_responses[atom] = isotropic
        anisotropic_responses[atom] = anisotropic

    for atom, amounts in  all_averages.items():
        isotropic: dict = {}
        anisotropic: dict = {}
        for label, averages in amounts.items():
            isotropic[label] = (averages[a] + averages[b] + averages[c])/3.0
            anisotropic[label] = sig_x*averages[a] + sig_y*averages[b] + sig_z*averages[c]
        isotropic_averages[atom] = isotropic
        anisotropic_averages[atom] = anisotropic

    return isotropic_responses, isotropic_averages, anisotropic_responses, anisotropic_averages

############## Print

def print_lresc_brief(io: scratch = None,
                      isotropic_responses: dict = None,
                      isotropic_averages: dict = None,
                      anisotropic_responses: dict = None,
                      anisotropic_averages: dict = None
                    ):
    """
    Brief information about LRESC's results

    Args:
    ----
        io (object:scratch): Driver to driver the output and binary files
        isotropic_responses (dict): Isotropic values of different lresc's responses
        isotropic_averages (dict): Isotropic values of different lresc's averages
        anisotropic_responses (dict): Anisotropic values of different lresc's responses
        anisotropic_averages (dict): Anisotropic values of different lresc's averages
    """
    # Brief results
    # -- isotropic
    paranr: float = isotropic_responses["paranr"]
    dianr: float = isotropic_averages["dianr"]
    lresc_para: float = sum([correction for label, correction in isotropic_responses.items()
                        if label not in ["paranr"]])
    lresc_dia: float = (sum([correction for label, correction in isotropic_averages.items()
                    if label not in ["dianr"]]))
    ligand_correction: float = sum([correction for label, correction in isotropic_responses.items()
                        if label in ["lpsokin", "lkinpso"]])
    core_correction: float = (sum([correction for label, correction in isotropic_averages.items()
                        if label not in ["dianr"]])
                        + sum([correction for label, correction in isotropic_responses.items()
                    if label not in ["lpsokin", "lkinpso"]]))
    # -- anisotropic
    ani_paranr: float = anisotropic_responses["paranr"]
    ani_dianr: float = anisotropic_averages["dianr"]
    ani_lresc_para = sum([correction for label, correction in anisotropic_responses.items()
                        if label not in ["paranr"]])
    ani_lresc_dia = (sum([correction for label, correction in anisotropic_averages.items()
                    if label not in ["dianr"]]))
    ani_ligand_correction = sum([correction for label, correction in anisotropic_responses.items()
                        if label in ["lpsokin", "lkinpso"]])
    ani_core_correction = (sum([correction for label, correction in anisotropic_averages.items()
                        if label not in ["dianr"]])
                        + sum([correction for label, correction in anisotropic_responses.items()
                    if label not in ["lpsokin", "lkinpso"]]))

    io.write_output("    " + "NR".center(21) + "Corrections".center(21) + "Total".center(31))
    io.write_output("    " + "Para".center(10) + " " + "Dia".center(10) + " "
            "Para".center(10) + " " + "Dia".center(10) + " "
            "NR".center(10) + " " + "Corrections".center(10) + " "
            "LRESC".center(10))
    io.write_output("-"*80)
    # Iso
    io.write_output("iso " + f"{paranr:.4f}".center(10) + " " + f"{dianr:.4f}".center(10) + " "
    f"{lresc_para:.4f}".center(10) + " " + f"{lresc_dia:.4f}".center(10) + " "
    f"{paranr+dianr:.4f}".center(10) + " " + f"{lresc_para + lresc_dia:.4f}".center(10)
    + " " + f"{paranr+dianr+lresc_para+lresc_dia:.4f}".center(10))
    # Ani
    io.write_output("ani " + f"{ani_paranr:.4f}".center(10) + " " + f"{ani_dianr:.4f}".center(10) + " "
    f"{ani_lresc_para:.4f}".center(10) + " " + f"{ani_lresc_dia:.4f}".center(10) + " "
    f"{ani_paranr+ani_dianr:.4f}".center(10) + " " + f"{ani_lresc_para + ani_lresc_dia:.4f}".center(10)
    + " " + f"{ani_paranr+ani_dianr+ani_lresc_para+ani_lresc_dia:.4f}".center(10))

    io.write_output(" "*25 + "-"*21)
    io.write_output("    " + " "*21 + "Ligand".center(10) + " " + "Core".center(10))
    io.write_output(" "*25 + "-"*21)
    io.write_output("iso " + " "*21 + f"{ligand_correction:.4f}".center(10) + " " + f"{core_correction:.4f}".center(10))
    io.write_output("ani " + " "*21 + f"{ani_ligand_correction:.4f}".center(10) + " " + f"{ani_core_correction:.4f}".center(10))


def print_lresc_values( io: scratch = None,
                        isotropic_responses: dict = None, 
                        isotropic_averages: dict = None,
                        anisotropic_responses: dict = None, 
                        anisotropic_averages: dict = None,
                        atom_label: list = None,
                        verbose: int = 0):
    """
    Driver print isotropic and anisotropic results

    Args:
    ----
        io (object:scratch): Driver to driver the output and binary files
        isotropic_responses (dict): Isotropic values of different lresc's responses
                                    of all atoms
        isotropic_averages (dict): Isotropic values of different lresc's averages
                                    of all atoms
        anisotropic_responses (dict): Anisotropic values of different lresc's responses
                                    of all atoms
        anisotropic_averages (dict): Anisotropic values of different lresc's averages
                                    of all atoms
        atom_label (list): Atoms' number
    """

    for atom in isotropic_responses.keys():

        if atom > 0:
            io.write_output(":"*89)
        io.write_output(f"@@@@ Atom: {str(atom_label[atom]) + str(atom + 1):} @@@@".center(89))

        if verbose > 10:

            if "paranr" in isotropic_responses[atom].keys():
                io.write_output("---> Non-Relativistic Paramagnetic <---".center(76))
                io.write_output(information=[lresc_label["paranr"]], type = 4,
                                values=[isotropic_responses[atom]["paranr"],
                                        anisotropic_responses[atom]["paranr"]])
            # Paramagnetic corrections
            names: list = [lresc_label[name] for name in isotropic_responses[atom].keys()
                            if name != "paranr" and name != "a2mv" and name != "a2dw"]
            values_iso: list = [value for name, value in isotropic_responses[atom].items()
                                if name != "paranr" and name != "a2mv" and name != "a2dw"]
            values_ani: list = [value for name, value in anisotropic_responses[atom].items()
                                if name != "paranr" and name != "a2mv" and name != "a2dw"]
            if len(names) > 0:
                io.write_output("---> Paramagnetic Corrections <---".center(76))
                io.write_output(information=names, values=values_iso + values_ani, type = 4)

            if "dianr" in isotropic_averages[atom].keys():
                io.write_output("---> Non-Relativistic Diamagnetic <---".center(76))
                io.write_output(information=[lresc_label["dianr"]], type = 4,
                                values=[isotropic_averages[atom]["dianr"],
                                        anisotropic_averages[atom]["dianr"]])
            # Diamagnetic corrections
            names: list = [lresc_label[name] for name in isotropic_averages[atom].keys()
                            if name != "dianr"]
            names += [lresc_label[name] for name in isotropic_responses[atom].keys()
                        if name == "a2mv" or name == "a2dw"]
            values_iso: list = [value for name, value in isotropic_averages[atom].items()
                                if name != "dianr"]
            values_iso += [value for name, value in isotropic_responses[atom].items()
                            if name == "a2mv" or name == "a2dw"]
            values_ani: list = [value for name, value in anisotropic_averages[atom].items()
                                if name != "dianr"]
            values_ani += [value for name, value in anisotropic_responses[atom].items()
                            if name == "a2mv" or name == "a2dw"]
            if len(names) > 0:
                io.write_output("---> Diamagnetic Corrections <---".center(76))
                io.write_output(information=names, values=values_iso + values_ani, type = 4)

        print_lresc_brief(
                        io = io,
                        isotropic_responses = isotropic_responses[atom],
                        isotropic_averages = isotropic_averages[atom],
                        anisotropic_responses = anisotropic_responses[atom],
                        anisotropic_averages = anisotropic_averages[atom]
                        )

def print_lresc_tensor(io: scratch = None,
                        responses_tensor: dict = None, 
                        averages_tensor: list = None,
                        isotropic_responses: dict = None, 
                        isotropic_averages: dict = None,
                        anisotropic_responses: dict = None, 
                        anisotropic_averages: dict = None,
                        atom_label: list = None):
    """
    Driver the tensor print

    Args:
    ----
        io (object:scratch): Driver to driver the output and binary files
        responses_tensor (dict): Tensor of different lresc's responses
        averages_tensor (list): Tensor of different lresc's averages
        isotropic_responses (dict): Isotropic values of different lresc's responses
                                    of all atoms
        isotropic_averages (dict): Isotropic values of different lresc's averages
                                    of all atoms
        anisotropic_responses (dict): Anisotropic values of different lresc's responses
                                    of all atoms
        anisotropic_averages (dict): Anisotropic values of different lresc's averages
                                    of all atoms
        atom_label (list): Atoms' number
    """
    for atom in responses_tensor.keys():

        if atom > 0:
            io.write_output(":"*89)
        io.write_output(f"@@@@ Atom: {str(atom_label[atom]) + str(atom + 1):} @@@@".center(101))

        if "paranr" in responses_tensor[atom].keys():
            io.write_output("---> Non-Relativistic Paramagnetic <---".center(101))
            io.write_output(information=[lresc_label["paranr"]], type = 4, box_type = 1,
                            values=[responses_tensor[atom]["paranr"]])

        # Paramagnetic corrections
        names: list = [lresc_label[name] for name in responses_tensor[atom].keys()
                        if name != "paranr" and name != "a2mv" and name != "a2dw"]
        values: list = [value for name, value in responses_tensor[atom].items()
                        if name != "paranr" and name != "a2mv" and name  != "a2dw"]
        if len(names) > 0:
            io.write_output("---> Paramagnetic Corrections <---".center(101))
            io.write_output(information=names, values=values, type=4, box_type=1)

        if "dianr" in averages_tensor[atom].keys():
            io.write_output("---> Non-Relativistic Diamagnetic <---".center(101))
            io.write_output(information=[lresc_label["dianr"]], type=4, box_type = 1,
                            values=[averages_tensor[atom]["dianr"]])

        # Diamagnetic corrections
        names: list = [lresc_label[name] for name in averages_tensor[atom].keys()
                        if name != "dianr"]
        names += [lresc_label[name] for name in responses_tensor[atom].keys()
                    if name == "a2mv" or name  == "a2dw"]
        values: list = [value for name, value in averages_tensor[atom].items()
                        if name != "dianr"]
        values += [value for name, value in responses_tensor[atom].items()
                    if name == "a2mv" or name  == "a2dw"]
        if len(names) > 0:
            io.write_output("---> Diamagnetic Corrections <---".center(101))
            io.write_output(information=names, values=values, type=4, box_type=1)

        print_lresc_brief(
                        io = io,
                        isotropic_responses = isotropic_responses[atom],
                        isotropic_averages = isotropic_averages[atom],
                        anisotropic_responses = anisotropic_responses[atom],
                        anisotropic_averages = anisotropic_averages[atom]
                        )
