from libp import *

def correction_to_calculate(ppesc_amounts: list = None):
    """
    Turn on or off ppesc_amounts corrections
    Args:
        paramagneitc (list): Paramagnetic amount to calculate
    """
    triplet_lineal_amount: bool = False
    singlet_lineal_amount: bool = False
    para_nr: bool = False
    dia_nr: bool = False
    diaavs: bool = False
    if not ppesc_amounts:
        return {
        "paranr": True,
        # Lineal -- Triplet
        "fclap": True, "sddxx": True, "sddyy": True, "sddzz": True,
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
        "sdkinyy": "sddyy" in ppesc_amounts,"sddzz": "sddzz" in ppesc_amounts,
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

def get_average(wf: wave_function = None, atom: int = None, driver_time: drv_time = None,
                gaugeo: list = None, type_correction: dict = None, corrections: dict = None,
                ppesc_consts: dict = None, verbose: int = 0, verbose_average: int = -1,
                ):
    """_summary_

    Args:
    ----
        wf (wave_function): _description_. Defaults to None.
        atom (int): _description_. Defaults to None.
        driver_time (drv_time): _description_. Defaults to None.
        gaugeo (list): _description_. Defaults to None.
        type_correction (dict): _description_. Defaults to None.
        corrections (dict): _description_. Defaults to None.
        verbose (int): _description_. Defaults to 0.
    """
    average_values: dict = {}
    av: average = average(wf)
    for name, amount in type_correction.items():
        if corrections[name]:
            value: dict = {}
            for label, operator in amount.items():
                value[label] = (ppesc_consts[name]*
                    list(av.calculate_average(
                    property = operator(atom),
                    gaugeo = gaugeo,
                    verbose = verbose_average
                    ).values())[0])
            average_values[name] = value
    return average_values

def get_gpv(wf: wave_function = None, atom: int = None, driver_time: drv_time = None,
            gaugeo: list = None, type_correction: dict = None, corrections: dict = None,
            quadratic: bool = False, verbose: int = 0, verbose_int: int = -1
            ):
    """
    Get all gradient properties and so on, for lresc calculation

    Args:
    ----
        gaugeo (list): Gauge
        atom (int): Atom index
        type_correction (dict): Type of correction
        corrections (dict): Dictionary that activate of corrections to calculate
        verbose_responses (int): Print level for response
    """
    gpvs: dict = {}
    for name, amount in type_correction.items():
        if corrections[name]:
            gpv: dict = {}
            mo_occupied: dict = {}
            mo_virtual: dict = {}
            for label, response in amount.items():
                mo_occupied[label], mo_virtual[label], gpv[label] = drv_gradient_property_vector(wf = wf,
                                        properties = response(atom),
                                        driver_time = driver_time, gaugeo = gaugeo,
                                        average = quadratic,
                                        verbose = verbose, verbose_integrals = verbose_int)
            gpvs[name] = gpv

    return gpvs

def get_lineal_response(responses: dict = None, type_correction: dict = None,
                        ppesc_consts: dict = None, rotations: int = None, corrections: dict = None,
                        rotation_i_a: int = None, rotation_j_b: int = None,
                        gpvs: dict = None, pp_left: float = None, pp_right: float = None):

    values: dict = {}
    for name, amount in type_correction.items():
        if corrections[name]:
            value: dict = {}
            for label in amount.keys():
                # it's multiplicated by minus the value like dalton do in the responses
                value[label] = 0.5*ppesc_consts[name]*(
                                    list(gpvs[name][label].values())[0][rotation_i_a]*
                                    pp_left*
                                    list(gpvs[name][label].values())[1][rotation_j_b+rotations]
                                +
                                    list(gpvs[name][label].values())[1][rotation_j_b]*
                                    pp_right*
                                    list(gpvs[name][label].values())[0][rotation_i_a+rotations]
                                )
                if name in responses.keys():
                    if label in responses[name].keys():
                        value[label] += responses[name][label]
            values[name] = value

    return values

def run_shielding_lresc(wf: wave_function = None, ppesc_amounts: list = None,
                                ppesc_consts: dict = None, atom: list = None,
                                principal_propagator_approximation: str = None,
                                driver_time: drv_time = None, verbose: int = 0,
                                verbose_integrals: int = -1):
    """
    Calculate all ppesc_amounts values
    Args:
        wf (wave_function): Wave function object
        paramagneitc (list): Paramagnetic amount to calculate
        ppesc_consts (dict): Dictionary with constants for differente amounts
        atom (list): Atom index to calculate the shielding
        principal_propagator_approximation (str): Approximation
        driver_time (drv_time): Object to manage the time calculation
        verbose (int): Print level
        verbose_integrals (int): Print level for integral module
    """
    start: float = time()

    corrections, responses_amounts, averages = correction_to_calculate(ppesc_amounts)
    label_amounts: dict = {"para_nr": paramagnetic_nr,
                    "triplet_lineal_amount": triplet_lineal_responses,
                    "singlet_lineal_amount": singlet_lineal_responses,
                    "dia_nr": diamagnetic_nr, "dia_avs": dia_averages}

    # parameters
    nocc: int = wf.mo_occ
    nvir: int = wf.mo_virt
    moe: list = wf.mo_energies

    # exchange and coulomb integrals
    if responses_amounts.values():
        coulomb_integrals, exchange_integrals = get_coulomb_exchange_integrals(wf,
                                            time_object = driver_time,
                                            verbose = verbose,
                                            verbose_int = verbose_integrals)
    # end exchange and coulomb

    # build principal propagator
    triplet: bool = responses_amounts["triplet_lineal_amount"]
    singlet: bool = (responses_amounts["para_nr"] or responses_amounts["singlet_lineal_amount"])
    if triplet or singlet:
        lineal_pp = drv_principal_propagator(driver_time = driver_time, moe = moe,
                                            n_mo_occ = nocc, n_mo_virt = nvir,
                                            coulomb = coulomb_integrals, exchange = exchange_integrals,
                                            multiplicity_pp = {"singlet": singlet,
                                                                "triplet": triplet},
                                            quadratic = False,
                                            tp_inv = 0, verbose = verbose)
    # end principal propagator

    gpvs: dict = {}
    occ_op_occ: dict = {}
    vir_op_vir: dict = {}
    all_averages: dict = {}
    for a in atom:
        gaugeo: list = wf.coordinates[a]

        atom_av: dict = {}
        if True in averages.values():
            for label, activate in averages.items():
                if activate:
                    temp_av = get_average(
                        wf = wf, atom = a, driver_time = driver_time, gaugeo = gaugeo,
                        type_correction = label_amounts[label], corrections = corrections,
                        ppesc_consts = ppesc_consts, verbose = verbose
                    )
                if activate: atom_av.update(temp_av)
        all_averages[a] = atom_av

        atom_gpvs: dict = {}
        if singlet or triplet:
            for label, activate in responses_amounts.items():
                temp_gpv = get_gpv(
                        wf = wf, atom = a, driver_time = driver_time, gaugeo = gaugeo,
                        type_correction = label_amounts[label], corrections = corrections,
                        quadratic = "False", verbose = verbose, verbose_int = verbose_integrals
                    )
                if activate:
                    atom_gpvs.update(temp_gpv)
            gpvs[a] = atom_gpvs

    # Responses
    all_responses: dict = {}
    rotations: int = nocc*nvir
    for i in range(nocc):
        for a in range(nvir):
            s = a + nocc
            for j in range(nocc):
                for b in range(nvir):
                    t = b + nocc

                    if singlet or triplet:
                        for at in atom:
                            atom_responses: dict = {}
                            for label, activate in responses_amounts.items():
                                if activate and ("lineal" in label or "para_nr" == label):
                                    if "triplet" in label:
                                        pp_left: float = lineal_pp["triplet"][a+i*nvir,b+j*nvir]
                                        pp_right: float = lineal_pp["triplet"][b+j*nvir,a+i*nvir]
                                    else:
                                        pp_left: float = lineal_pp["singlet"][a+i*nvir,b+j*nvir]
                                        pp_right: float = lineal_pp["singlet"][b+j*nvir,a+i*nvir]
                                    if at not in all_responses.keys(): path_responses: dict = {}
                                    else: path_responses: dict = all_responses[at]
                                    temp_value: dict = {}
                                    temp_value = get_lineal_response(
                                        responses = path_responses, type_correction = label_amounts[label],
                                        rotations = rotations, corrections = corrections,
                                        ppesc_consts = ppesc_consts, rotation_i_a = a+i*nvir,
                                        rotation_j_b = b+j*nvir, gpvs = gpvs[at],
                                        pp_left = pp_left, pp_right = pp_right
                                    )
                                if activate: atom_responses.update(temp_value)
                            all_responses[at] = atom_responses

    if verbose > 10:
        driver_time.add_name_delta_time(name = "PPESC Amounts Calculations", delta_time = (time() - start))

    return all_responses, all_averages

def get_shielding_isotropic(all_responses: dict = None, all_averages: dict = None):

    isotropic_responses: dict = {}
    isotropic_averages: dict = {}

    if all_responses:
        isotropic_atom: dict = None
        for atom, dict_corrections in  all_responses.items():
            isotropic: dict = {}
            for correction, components_values in dict_corrections.items():
                isotropic[correction] = sum(components_values.values())/3.0
            #
            if not isotropic_atom: isotropic_atom = isotropic
            else: isotropic_atom.update(isotropic)
            #
            isotropic_responses[atom] = isotropic_atom
            if "sddxx" in isotropic_responses[atom].keys():
                isotropic_responses[atom]["sdlap"] = (isotropic_responses[atom]["sddxx"] +
                                isotropic_responses[atom]["sddyy"] +isotropic_responses[atom]["sddzz"])

    if all_averages:
        isotropic_atom: dict = None
        for atom, dict_corrections in  all_averages.items():
            isotropic: dict = {}
            for correction, components_values in dict_corrections.items():
                if "fc" != correction:
                    isotropic[correction] = sum(components_values.values())/3.0
                else:
                    isotropic["fc"] = sum(components_values.values())
            #
            if not isotropic_atom: isotropic_atom = isotropic
            else: isotropic_atom.update(isotropic)
            #
            isotropic_averages[atom] = isotropic_atom

    return isotropic_responses, isotropic_averages

def get_shielding_anisotropic(all_responses: dict = None, all_averages: dict = None):

    anisotropic_responses: dict = {}
    anisotropic_averages: dict = {}

    # Anisotropic with z
    z_component: list = [
                            "lzpso3", "fcdzz", "angmomzpsoke3", "psozozke3",
                            "nstcgo3z", "psooz3z", "dnske3z", "sd3z", "pnstcgop3z"
                        ]
    z_sign: float = -1.0

    if all_responses:
        anisotropic_atom: dict = None
        for atom, dict_corrections in  all_responses.items():
            anisotropic: dict = {}
            for correction, components_values in dict_corrections.items():
                if len(components_values) == 3 and correction not in ["sddxx", "sddyy", "sddzz"]:
                    anisotropic[correction] = sum([value if name not in z_component else -value
                                                for name, value in components_values.items()])
                else:
                    anisotropic[correction] = sum(components_values.values())

            #
            if not anisotropic_atom: anisotropic_atom = anisotropic
            else: anisotropic_atom.update(anisotropic)
            #
            anisotropic_responses[atom] = anisotropic_atom
            if "sddxx" in anisotropic_responses[atom].keys():
                anisotropic_responses[atom]["sdlap"] = (anisotropic_responses[atom]["sddxx"] +
                                anisotropic_responses[atom]["sddyy"] + z_sign*anisotropic_responses[atom]["sddzz"])

    if all_averages:
        anisotropic_atom: dict = None
        for atom, dict_corrections in  all_averages.items():
            anisotropic: dict = {}
            for correction, components_values in dict_corrections.items():
                if correction != "fc":
                    anisotropic[correction] = sum([value if name not in z_component else -value
                                                for name, value in components_values.items()])
                else:
                    anisotropic["fc"] = sum(components_values.values())
            #
            if not anisotropic_atom: anisotropic_atom = anisotropic
            else: anisotropic_atom.update(anisotropic)
            #
            anisotropic_averages[atom] = anisotropic_atom

    return anisotropic_responses, anisotropic_averages


############## Print

def print_sinlget_lineal(responses: dict = None, ani_responses: dict = None):
    if "lpsokin" in responses.keys():
        lpsokin: float = responses["lpsokin"]
        ani_lpsokin: float = ani_responses["lpsokin"]
    else:
        lpsokin: str = "No Calc."
        ani_lpsokin: str = "No Calc."
    if "lkinpso" in responses.keys():
        lkinpso: float = responses["lkinpso"]
        ani_lkinpso: float = ani_responses["lkinpso"]
    else:
        lkinpso: str = "No Calc."
        ani_lkinpso: str = "No Calc."

    print(("<<L; {PSO, K}>>".center(20) + "   " + "<<PSO; {L, K}>>".center(20)).center(89))
    print(("|" + "¯"*18 + "|" + "   " + "|" + "¯"*18 + "|").center(89))
    print(("|" + f"{lpsokin:.6f}".center(18) + "|" + "   " + "|" + f"{lkinpso:.6f}".center(18) + "|").center(89))
    print(("|" + f"{ani_lpsokin:.6f}".center(18) + "|" + "   " + "|" + f"{ani_lkinpso:.6f}".center(18) + "|").center(89))
    print(("¯"*18 + "     " + "¯"*18).center(89))

def print_triplet_lineal(responses: dict = None, ani_responses: dict = None):
    if "fclap" in responses.keys():
        fckin: float = responses["fclap"]
        ani_fckin: float = ani_responses["fclap"]
    else:
        fckin: str = "No Calc."
        ani_fckin: str = "No Calc."
    if "sdlap" in responses.keys():
        sdkin: float = responses["sdlap"]
        ani_sdkin: float = ani_responses["sdlap"]
    else:
        sdkin: str = "No Calc."
        ani_sdkin: str = "No Calc."
    print(("<<FC; K>>".center(20) + "   " + "<<SD; K>>".center(20)).center(89))
    print(("|" + "¯"*18 + "|" + "   " + "|" + "¯"*18 + "|").center(89))
    print(("|" + f"{fckin:.6f}".center(18) + "|" + "   " + "|" + f"{sdkin:.6f}".center(18) + "|").center(89))
    print(("|" + f"{ani_fckin:.6f}".center(18) + "|" + "   " + "|" + f"{ani_sdkin:.6f}".center(18) + "|").center(89))
    print(("¯"*18 + "     " + "¯"*18).center(89))


def print_dia_averages(averages: dict = None, ani_averages: dict = None):
    if "fc" in averages.keys():
        fc: float = averages["fc"]
        ani_fc: float = ani_averages["fc"]
    else:
        fc: str = "No Calc."
        ani_fc: str = "No Calc."
    if "psooz" in averages.keys():
        psooz: float = averages["psooz"]
        ani_psooz: float = averages["psooz"]
    else:
        psooz: str = "No Calc."
        ani_psooz: str = "No Calc."
    if "dnske" in averages.keys():
        dnske: float = averages["dnske"]
        ani_dnske: float = averages["dnske"]
    else:
        dnske: str = "No Calc."
        ani_dnske: str = "No Calc."
    if "sd" in averages.keys():
        sd: float = averages["sd"]
        ani_sd: float = averages["sd"]
    else:
        sd: str = "No Calc."
        ani_sd: str = "No Calc."
    if "pnstcgop" in averages.keys():
        pnstcgop: float = averages["pnstcgop"]
        ani_pnstcgop: float = averages["pnstcgop"]
    else:
        pnstcgop: str = "No Calc."
        ani_pnstcgop: str = "No Calc."

    print(("<FC>".center(20) + "   " + "<PSOOZ>".center(20) + "   " + "<DNSKE>".center(20)).center(89))
    print(("|" + "¯"*18 + "|" + "   " + "|" + "¯"*18 + "|" + "   " + "|" + "¯"*18 + "|").center(89))
    print(("|" + f"{fc:.6f}".center(18) + "|" + "   " + "|" + f"{psooz:.6f}".center(18) + "|"
            + "   " + "|" + f"{dnske:.6f}".center(18) + "|").center(89))
    print(("|" + f"{ani_fc:.6f}".center(18) + "|" + "   " + "|" + f"{ani_psooz:.6f}".center(18) + "|"
            + "   " + "|" + f"{ani_dnske:.6f}".center(18) + "|").center(89))
    print(("¯"*18 + "     " + "¯"*18  + "     " + "¯"*18).center(89))
    print()
    print(("<SD>".center(20) + "   " + "<PNSTCGOP>".center(20)).center(89))
    print(("|" + "¯"*18 + "|" + "   " + "|" + "¯"*18 + "|").center(89))
    print(("|" + f"{sd:.6f}".center(18) + "|" + "   " + "|" + f"{pnstcgop:.6f}".center(18) + "|").center(89))
    print(("|" + f"{ani_sd:.6f}".center(18) + "|" + "   " + "|" + f"{ani_pnstcgop:.6f}".center(18) + "|").center(89))
    print(("¯"*18 + "     " + "¯"*18).center(89))


def print_lresc_values(isotropic_responses: dict = None, isotropic_averages: dict = None,
                        anisotropic_responses: dict = None, anisotropic_averages: dict = None,
                        atom_label: list = None):
    """
    Driver print isotropic and anisotropic results

    Args:
        isotropic_responses (dict): 
        isotropic_averages (dict): 
        anisotropic_averages (dict): 
        anisotropic_responses (dict): 
        atom_label (list): 
    """
    print()
    print(("*"*27).center(70))
    print(f"* LRESC Shielding Results *".center(70))
    print(("*"*27).center(70))

    print(f"Isotropic Values [ppm]".center(89), "\n", "and".center(89), "\n", "Anisotropic Values [ppm]".center(89))
    print(("-"*40).center(89))

    for count, atom in enumerate(isotropic_responses.keys()):

        if count > 0:
            print("\n",":"*89)
        print("\n\n",f"@@@@ Atom: {str(atom_label[atom]) + str(atom + 1):} @@@@".center(89),"\n\n")

        if "paranr" in isotropic_responses[atom].keys():
            paranr: float = isotropic_responses[atom]["paranr"]
            ani_paranr: float = anisotropic_responses[atom]["paranr"]
        else:
            paranr: str = "No Calc."
            ani_paranr: str = "No Calc."
        if "paranr" in isotropic_responses[atom].keys():
            print((f"Para".center(20)).center(89))
            print(("|" + "¯"*18 + "|").center(89))
            print(("|" + f"{paranr:.6f}".center(18) + "|").center(89))
            print(("|" + f"{ani_paranr:.6f}".center(18) + "|").center(89))
            print(("¯"*18).center(89))

    # Paramagnetic corrections
        print("\n","---> Paramagnetic Corrections <---".center(89),"\n")
        print_sinlget_lineal(isotropic_responses[atom],anisotropic_responses[atom])
        print_triplet_lineal(isotropic_responses[atom],anisotropic_responses[atom])

        if "dianr" in isotropic_averages[atom].keys():
            dianr: float = isotropic_averages[atom]["dianr"]
            ani_dianr: float = anisotropic_averages[atom]["dianr"]
        else:
            dianr: str = "No Calc."
            ani_dianr: str = "No Calc."
        if "dianr" in isotropic_averages[atom].keys():
            print((f"Dia".center(20)).center(89))
            print(("|" + "¯"*18 + "|").center(89))
            print(("|" + f"{dianr:.6f}".center(18) + "|").center(89))
            print(("|" + f"{ani_dianr:.6f}".center(18) + "|").center(89))
            print(("¯"*18).center(89))

    # Diamagnetic corrections
        print("\n","---> Diamagnetic Corrections <---".center(89),"\n")
        print_dia_averages(isotropic_averages[atom], anisotropic_averages[atom])

    # Brief results
        # -- isotropic
        lresc_para = sum([correction for label, correction in isotropic_responses[atom].items()
                            if label not in ["paranr"]])
        lresc_dia = (sum([correction for label, correction in isotropic_averages[atom].items()
                        if label not in ["dianr"]]))
        ligand_correction = sum([correction for label, correction in isotropic_responses[atom].items()
                            if label in ["lpsokin", "lkinpso"]])
        core_correction = (sum([correction for label, correction in isotropic_averages[atom].items()
                            if label not in ["dianr"]])
                            + sum([correction for label, correction in isotropic_responses[atom].items()
                        if label not in ["lpsokin", "lkinpso"]]))
        # -- anisotropic
        ani_lresc_para = sum([correction for label, correction in anisotropic_responses[atom].items()
                            if label not in ["paranr"]])
        ani_lresc_dia = (sum([correction for label, correction in anisotropic_averages[atom].items()
                        if label not in ["dianr"]]))
        ani_ligand_correction = sum([correction for label, correction in anisotropic_responses[atom].items()
                            if label in ["lpsokin", "lkinpso"]])
        ani_core_correction = (sum([correction for label, correction in anisotropic_averages[atom].items()
                            if label not in ["dianr"]])
                            + sum([correction for label, correction in anisotropic_responses[atom].items()
                        if label not in ["lpsokin", "lkinpso"]]))

        print("    " + "NR".center(21) + "LRESC".center(21))
        print("    " + "Para".center(10) + " " + "Dia".center(10) + " "
                "Para".center(10) + " " + "Dia".center(10) + " "
                "NR".center(10) + " " + "LRESC".center(10) + " "
                "Total".center(10))
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
