from urllib import response
from libp import *

def correction_to_calculate(ppesc_amounts: list = None, tensor: bool = False):
    """
    Turn on or off ppesc_amounts corrections
    Args:
        paramagneitc (list): Paramagnetic amount to calculate
        tensor (bool): Activate calculation of all tensor components
    """
    triplet_lineal_amount: bool = False
    singlet_lineal_amount: bool = False
    para_nr: bool = False
    dia_nr: bool = False
    diaavs: bool = False

    sddxy: bool = False
    sddxz: bool = False
    sddyx: bool = False
    sddyz: bool = False
    sddzx: bool = False
    sddzy: bool = False
    if tensor:
        sddxy: bool = True
        sddxz: bool = True
        sddyx: bool = True
        sddyz: bool = True
        sddzx: bool = True
        sddzy: bool = True
    if not ppesc_amounts:
        return {
        "paranr": True,
        # Lineal -- Triplet
        "fclap": True, "sddxx": True, "sddyy": True, "sddzz": True,
        "sddxy": sddxy, "sddxz": sddxz, "sddyx": sddyx, "sddyz": sddyz,
        "sddzx": sddzx, "sddzy": sddzy,
        # Lineal -- Singlet
        "lpsokin": True, "lkinpso": True,
        # Diamagnetic
        # NR
        "dianr": True,
        # Averages
        "fc": True, "sd": True, "psooz": True, "dnske": True, "pnstcgop": True
        }, {"para_nr": True, "triplet_lineal_amount": True, "singlet_lineal_amount": True,
        },{"dia_nr": True, "dia_avs": True}
    else:
        if "paranr" in ppesc_amounts or "nr" in ppesc_amounts: para_nr = True
        if "dianr" in ppesc_amounts or "nr" in ppesc_amounts: dia_nr = True
        if "diaavs" in ppesc_amounts or "avs" in ppesc_amounts: diaavs = True
        if ("fclap" in ppesc_amounts or "sdlap" in ppesc_amounts): triplet_lineal_amount = True
        if "lpsokin" in ppesc_amounts or "lkinpso" in ppesc_amounts: singlet_lineal_amount = True
        return {
        "paranr": "paranr" in ppesc_amounts,
        # Lineal -- Triplet
        "fclap": "fclap" in ppesc_amounts, "sddxx": "sddxx" in ppesc_amounts,
        "sddyy": "sddyy" in ppesc_amounts,"sddzz": "sddzz" in ppesc_amounts,
        "sddxy": "sddxy" in ppesc_amounts, "sddxz": "sddxz" in ppesc_amounts,
        "sddyx": "sddyx" in ppesc_amounts, "sddyz": "sddyz" in ppesc_amounts,
        "sddzx": "sddzx" in ppesc_amounts, "sddzy": "sddzy" in ppesc_amounts,
        # Lineal -- Singlet
        "lpsokin": "lpsokin" in ppesc_amounts, "lkinpso": "lkinpso" in ppesc_amounts,
        # Diamagnetic
        # NR
        "dianr": "dianr" in ppesc_amounts,
        # Averages
        "fc": "fc" in ppesc_amounts, "sd": "sd" in ppesc_amounts, "psooz": "psooz" in ppesc_amounts,
        "dnske": "dnske" in ppesc_amounts, "pnstcgop": "pnstcgop" in ppesc_amounts
        }, {"para_nr": para_nr, "triplet_lineal_amount": triplet_lineal_amount, "singlet_lineal_amount": singlet_lineal_amount,
        }, {"dia_nr": dia_nr, "dia_avs": diaavs}

def run_shielding(wf: wave_function = None, ppesc_amounts: list = None,
                    ppesc_consts: dict = None, atom: list = None,
                    principal_propagator_approximation: str = None,
                    driver_time: drv_time = None, verbose: int = 0,
                    tensor: bool = False,
                    scalar_correction: bool = None,
                    verbose_response: int = -1,
                    verbose_average: int = -1,
                    verbose_fock: int = 1):
    """
    Calculate all ppesc_amounts values

    The different amount involve in the PPESC methodology for Shielding
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
                                {"fclap": fclap, "sddxx": sddxx, "sddyy": sddyy, "sddzz": sddzz}
        singlet_lineal_amount: Correction type singlet response associated with
                                paramagnetic component
                                {"lpsokin": lpsokin, "lkinpso": lkinpso}
        dia_av: Correction type average associated with diamagnetic component
                {"fc": fc, "sd": sd, "psooz": psooz, "dnske": dnske, "pnstcgop": pnstcgop}
        }

    Also, there are another dictionary, which indicate if the amount will be calculated, it is
    called correction

    Args:
        wf (wave_function): Wave function object
        ppesc_amounts (list): Paramagnetic amount to calculate
        ppesc_consts (dict): Dictionary with constants for differente amounts
        atom (list): Atom index to calculate the shielding
        principal_propagator_approximation (str): Approximation
        tensor (bool): Activate calculation of all tensor components
        driver_time (drv_time): Object to manage the time calculation
        scalar_correction (bool): Activate Dw and Mv correction for energy
        verbose (int): Print level
        verbose_integrals (int): Print level for integral module
    """
    start: float = time()

    corrections, responses_amounts, averages = correction_to_calculate(ppesc_amounts, tensor)
    ## The values of next dictionary are in include ppesc_parameters.py
    label_amounts: dict = {"para_nr": paramagnetic_nr,
                    "triplet_lineal_amount": triplet_lineal_responses,
                    "singlet_lineal_amount": singlet_lineal_responses,
                    "dia_nr": diamagnetic_nr, "dia_avs": dia_averages}

    #Hidden Print of other object
    hidden_prints_other_object: object = HiddenPrints()
    hidden_prints_other_object.__enter__()
    # parameters
    if scalar_correction:
        if verbose_fock > 0: hidden_prints_other_object.__exit__(True,True,True)
        moe, at_2_in = fock().calculate_hf_moe(wf,
                                            at2in = True,
                                            relativity_correction=True,
                                            verbose=verbose_fock)

        if verbose_fock <= 0 and verbose_response > 0:
            hidden_prints_other_object.__exit__(True,True,True)
        else:
            hidden_prints_other_object.__enter__()

        lineal_response: response = response(wf = wf, moe = moe, at2in = at_2_in)
        if verbose_response < 0: hidden_prints_other_object.__exit__(True,True,True)
        driver_time.add_name_delta_time(name = "One-Body Mv and Dw Corrections to Energy", delta_time = (time() - start))
    else:
        moe: list = wf.mo_energies
        if verbose_response > 0: hidden_prints_other_object.__exit__(True,True,True)
        lineal_response: response = response(wf = wf, moe = moe)
        if verbose_response < 0: hidden_prints_other_object.__exit__(True,True,True)

    all_averages: dict = {}
    all_responses: dict = {}
    av: average = average(wf)

    delta_average: float = 0.0
    delta_response: float = 0.0
    for a in atom:
        gaugeo: list = wf.coordinates[a]

        # Avarage calculation
        start_average: float = time()
        if verbose_response <= 0 and verbose_average > 0:
            hidden_prints_other_object.__exit__(True,True,True)
        else:
            hidden_prints_other_object.__enter__()
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
                        temp_av_a: list = [(ppesc_consts[name]*
                                        list(av.calculate_average(
                                        property = operators[component](a),
                                        gaugeo = gaugeo,
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
        if verbose_average < 0: hidden_prints_other_object.__exit__(True,True,True)
        delta_average += time() - start_average
        # End Average calculation

        # Response calculation
        start_response: float = time()
        if verbose_average <= 0 and verbose_response > 0:
            hidden_prints_other_object.__exit__(True,True,True)
        else:
            hidden_prints_other_object.__enter__()
        atom_responses: dict = {}
        sdlap_components: list = ["sddxx", "sddxy", "sddxz", "sddyx", "sddyy", "sddyz", "sddzx", "sddzy", "sddzz"]
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
                        if label not in sdlap_components and label != "fclap":
                            tensor_components: int = 9
                        else:
                            tensor_components: int = 3
                            tensor_step: int = 1
                        temp_responses_a: list = [-ppesc_consts[label]
                                    *list(lineal_response.drv_reponse_calculation(
                                    gaugeo = gaugeo,
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
            if "fclap" in atom_responses.keys():
                temp: float = atom_responses["fclap"]
                atom_responses["fclap"] = [temp[0],0.0,0.0,0.0,temp[1],0.0,0.0,0.0,temp[2]]
        # Remove components of sdlap
        all_responses[a] =  {
                                name: atom_responses[name]
                                for name in name_order_responses
                            }
        delta_response += time() - start_response
        if verbose_response < 0: hidden_prints_other_object.__exit__(True,True,True)
        # End Response calculation

    if verbose > 10:
        driver_time.add_name_delta_time(name = "Averages Amounts Calculations", delta_time = delta_average)
        driver_time.add_name_delta_time(name = "Responses Amounts Calculations", delta_time = delta_response)
        driver_time.add_name_delta_time(name = "PPESC Amounts Calculations", delta_time = (time() - start))

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

def print_ppesc_brief(isotropic_responses: dict = None, isotropic_averages: dict = None,
                    anisotropic_responses: dict = None, anisotropic_averages: dict = None
                    ):
    """
    Brief information about LRESC's results

    Args:
    ----
        isotropic_responses (dict): Isotropic values of different ppesc's responses
        isotropic_averages (dict): Isotropic values of different ppesc's averages
        anisotropic_responses (dict): Anisotropic values of different ppesc's responses
        anisotropic_averages (dict): Anisotropic values of different ppesc's averages
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

    print("    " + "NR".center(21) + "Corrections".center(21) + "Total".center(31))
    print("    " + "Para".center(10) + " " + "Dia".center(10) + " "
            "Para".center(10) + " " + "Dia".center(10) + " "
            "NR".center(10) + " " + "Corrections".center(10) + " "
            "PPESC".center(10))
    print("-"*80)
    # Iso
    print("iso " + f"{paranr:.4f}".center(10) + " " + f"{dianr:.4f}".center(10) + " "
    f"{lresc_para:.4f}".center(10) + " " + f"{lresc_dia:.4f}".center(10) + " "
    f"{paranr+dianr:.4f}".center(10) + " " + f"{lresc_para + lresc_dia:.4f}".center(10)
    + " " + f"{paranr+dianr+lresc_para+lresc_dia:.4f}".center(10))
    # Ani
    print("ani " + f"{ani_paranr:.4f}".center(10) + " " + f"{ani_dianr:.4f}".center(10) + " "
    f"{ani_lresc_para:.4f}".center(10) + " " + f"{ani_lresc_dia:.4f}".center(10) + " "
    f"{ani_paranr+ani_dianr:.4f}".center(10) + " " + f"{ani_lresc_para + ani_lresc_dia:.4f}".center(10)
    + " " + f"{ani_paranr+ani_dianr+ani_lresc_para+ani_lresc_dia:.4f}".center(10))

    print(" "*25 + "-"*21)
    print("    " + " "*21 + "Ligand".center(10) + " " + "Core".center(10))
    print(" "*25 + "-"*21)
    print("iso " + " "*21 + f"{ligand_correction:.4f}".center(10) + " " + f"{core_correction:.4f}".center(10))
    print("ani " + " "*21 + f"{ani_ligand_correction:.4f}".center(10) + " " + f"{ani_core_correction:.4f}".center(10))

def print_ppesc_values(isotropic_responses: dict = None, isotropic_averages: dict = None,
                        anisotropic_responses: dict = None, anisotropic_averages: dict = None,
                        atom_label: list = None, verbose: int = 0):
    """
    Driver print isotropic and anisotropic results

    Args:
    ----
        isotropic_responses (dict): Isotropic values of different ppesc's responses
                                    of all atoms
        isotropic_averages (dict): Isotropic values of different ppesc's averages
                                    of all atoms
        anisotropic_responses (dict): Anisotropic values of different ppesc's responses
                                    of all atoms
        anisotropic_averages (dict): Anisotropic values of different ppesc's averages
                                    of all atoms
        atom_label (list): Atoms' number
    """

    for atom in isotropic_responses.keys():

        if atom > 0:
            print("\n",":"*89)
        print("\n\n",f"@@@@ Atom: {str(atom_label[atom]) + str(atom + 1):} @@@@".center(89),"\n\n")

        if verbose > 10:

            if "paranr" in isotropic_responses[atom].keys():
                print("\n","---> Non-Relativistic Paramagnetic <---".center(76),"\n")
                print_box(names=[ppesc_label["paranr"]],
                        values=[isotropic_responses[atom]["paranr"],
                        anisotropic_responses[atom]["paranr"]])

            # Paramagnetic corrections
            names: list = [ppesc_label[name] for name in isotropic_responses[atom].keys() if name != "paranr"]
            values_iso: list = [value for name, value in isotropic_responses[atom].items() if name != "paranr"]
            values_ani: list = [value for name, value in anisotropic_responses[atom].items() if name != "paranr"]
            if len(names) > 0:
                print("\n","---> Paramagnetic Corrections <---".center(76),"\n")
                print_box(names=names, values=values_iso + values_ani)

            if "dianr" in isotropic_averages[atom].keys():
                print("\n","---> Non-Relativistic Diamagnetic <---".center(76),"\n")
                print_box(names=[ppesc_label["dianr"]],
                        values=[isotropic_averages[atom]["dianr"],
                        anisotropic_averages[atom]["dianr"]])

            # Diamagnetic corrections
            names: list = [ppesc_label[name] for name in isotropic_averages[atom].keys() if name != "dianr"]
            values_iso: list = [value for name, value in isotropic_averages[atom].items() if name != "dianr"]
            values_ani: list = [value for name, value in anisotropic_averages[atom].items() if name != "dianr"]
            if len(names) > 0:
                print("\n","---> Diamagnetic Corrections <---".center(76),"\n")
                print_box(names=names, values=values_iso + values_ani)

        print_ppesc_brief(isotropic_responses = isotropic_responses[atom],
                        isotropic_averages = isotropic_averages[atom],
                        anisotropic_responses = anisotropic_responses[atom],
                        anisotropic_averages = anisotropic_averages[atom]
                        )

def print_ppesc_tensor(responses_tensor: dict = None, averages_tensor: list = None,
                        isotropic_responses: dict = None, isotropic_averages: dict = None,
                        anisotropic_responses: dict = None, anisotropic_averages: dict = None,
                        atom_label: list = None):
    """
    Driver the tensor print

    Args:
    ----
        responses_tensor (dict): Tensor of different ppesc's responses
        averages_tensor (list): Tensor of different ppesc's averages
        isotropic_responses (dict): Isotropic values of different ppesc's responses
                                    of all atoms
        isotropic_averages (dict): Isotropic values of different ppesc's averages
                                    of all atoms
        anisotropic_responses (dict): Anisotropic values of different ppesc's responses
                                    of all atoms
        anisotropic_averages (dict): Anisotropic values of different ppesc's averages
                                    of all atoms
        atom_label (list): Atoms' number
    """
    for atom in responses_tensor.keys():

        if atom > 0:
            print("\n",":"*89)
        print("\n\n",f"@@@@ Atom: {str(atom_label[atom]) + str(atom + 1):} @@@@".center(101),"\n\n")

        if "paranr" in responses_tensor[atom].keys():
            print("\n","---> Non-Relativistic Paramagnetic <---".center(101),"\n")
            print_tensor(names=[ppesc_label["paranr"]],
                    values=[responses_tensor[atom]["paranr"]])

        # Paramagnetic corrections
        names: list = [ppesc_label[name] for name in responses_tensor[atom].keys() if name != "paranr"]
        values: list = [value for name, value in responses_tensor[atom].items() if name != "paranr"]
        if len(names) > 0:
            print("\n","---> Paramagnetic Corrections <---".center(101),"\n")
            print_tensor(names=names, values=values)

        if "dianr" in averages_tensor[atom].keys():
            print("\n","---> Non-Relativistic Diamagnetic <---".center(101),"\n")
            print_tensor(names=[ppesc_label["dianr"]],
                    values=[averages_tensor[atom]["dianr"]])

        # Diamagnetic corrections
        names: list = [ppesc_label[name] for name in averages_tensor[atom].keys() if name != "dianr"]
        values: list = [value for name, value in averages_tensor[atom].items() if name != "dianr"]
        if len(names) > 0:
            print("\n","---> Diamagnetic Corrections <---".center(101),"\n")
            print_tensor(names=names, values=values)

        print_ppesc_brief(isotropic_responses = isotropic_responses[atom],
                        isotropic_averages = isotropic_averages[atom],
                        anisotropic_responses = anisotropic_responses[atom],
                        anisotropic_averages = anisotropic_averages[atom]
                        )