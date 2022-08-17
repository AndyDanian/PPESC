from optparse import Values
from libl import *

def correction_to_calculate(lresc_amounts: list = None):
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
    if not lresc_amounts:
        return {
        "paranr": True,
        # Lineal -- Triplet
        "fckin": True, "sdkinxx": True, "sdkinyy": True, "sdkinzz": True,
        "fcbso": True, "sdbsoxx": True, "sdbsoyy": True, "sdbsozz": True,
        # Lineal -- Singlet
        "lpsokin": True, "lkinpso": True,
        # Quadratic -- Triplet
        "lfcso" : True,
        "lsdsoxx": True, "lsdsoyy": True, "lsdsozz": True,
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
        if ("fckin" in lresc_amounts or "sdkin" in lresc_amounts or
            "fcbso" in lresc_amounts or "sdbso" in lresc_amounts): triplet_lineal_amount = True
        if "lpsokin" in lresc_amounts or "lkinpso" in lresc_amounts: singlet_lineal_amount = True
        if ("lfcso" in lresc_amounts or "lsdsoxx" in lresc_amounts or
            "lsdsoyy" in lresc_amounts or "lsdsozz" in lresc_amounts): triplet_quadratic_amount = True
        if ("lpsodw" in lresc_amounts or "lpsomv" in lresc_amounts): singlet_quadratic_amount = True
        if ("a2mv" in lresc_amounts or "a2dw" in lresc_amounts): dia_lineal = True
        return {
        "paranr": "paranr" in lresc_amounts,
        # Lineal -- Triplet
        "fckin": "fckin" in lresc_amounts, "sdkinxx": "sdkinxx" in lresc_amounts,
        "sdkinyy": "sdkinyy" in lresc_amounts,"sdkinzz": "sdkinzz" in lresc_amounts,
        "fcbso": "fcbso" in lresc_amounts, "sdbsoxx": "sdbsoxx" in lresc_amounts,
        "sdbsoyy": "sdbsoyy" in lresc_amounts, "sdbsozz": "sdbsozz" in lresc_amounts,
        # Lineal -- Singlet
        "lpsokin": "lpsokin" in lresc_amounts, "lkinpso": "lkinpso" in lresc_amounts,
        # Quadratic -- Triplet
        "lfcso" : "lfcso" in lresc_amounts,
        "lsdsoxx": "lsdsoxx" in lresc_amounts, "lsdsoyy": "lsdsoyy" in lresc_amounts, "lsdsozz": "lsdsozz" in lresc_amounts,
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

def get_average(wf: wave_function = None, atom: int = None, driver_time: drv_time = None,
                gaugeo: list = None, type_correction: dict = None, corrections: dict = None,
                lresc_consts: dict = None, verbose: int = 0, verbose_average: int = -1,
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
                value[label] = (lresc_consts[name]*
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
    mo_occupieds: dict = {}
    mo_virtuals: dict = {}
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
            mo_occupieds[name], mo_virtuals[name], gpvs[name] =  mo_occupied, mo_virtual, gpv

    return mo_occupieds, mo_virtuals, gpvs

def get_lineal_response(responses: dict = None, type_correction: dict = None,
                        lresc_consts: dict = None, rotations: int = None, corrections: dict = None,
                        rotation_i_a: int = None, rotation_j_b: int = None,
                        gpvs: dict = None, pp_left: float = None, pp_right: float = None):

    values: dict = {}
    for name, amount in type_correction.items():
        if corrections[name]:
            value: dict = {}
            for label in amount.keys():
                # it's multiplicated by minus the value like dalton do in the responses
                value[label] = -0.5*lresc_consts[name]*(
                                    list(gpvs[name][label].values())[0][rotation_i_a]*
                                    pp_left*
                                    list(gpvs[name][label].values())[1][rotation_j_b]
                                +
                                    list(gpvs[name][label].values())[1][rotation_j_b+rotations]*
                                    pp_right*
                                    list(gpvs[name][label].values())[0][rotation_i_a+rotations]
                                )
                if name in responses.keys():
                    if label in responses[name].keys():
                        value[label] += responses[name][label]
            values[name] = value

    return values

def get_quadratic_response(responses: dict = None, type_correction: dict = None,
                        lresc_consts: dict = None, rotations: int = None, corrections: dict = None,
                        rotation_i_a: int = None, rotation_j_b: int = None,
                        occ_occ: float = None, vir_vir: list = None, d_bc: bool = False,
                        occ_op_occ: dict = None, vir_op_vir: dict = None, d_ad: bool = False,
                        gpvs: dict = None, pp_op_a: list = None, pp_op_b: list = None):

    values: dict = {}
    for name, amount in type_correction.items():
        if corrections[name]:
            value: dict = {}
            for label in amount.keys():
                # it's multiplicated by minus the value like dalton do in the responses
                op_a: list = list(gpvs[name][label].values())[0]
                op_b: list = list(gpvs[name][label].values())[1]
                op_c: list = list(gpvs[name][label].values())[2]
                vir_a_occ: float = list(vir_op_vir[name][label].values())[0][vir_vir[0]]
                vir_1b_occ: float = list(vir_op_vir[name][label].values())[1][vir_vir[0]]
                vir_2b_occ: float = list(vir_op_vir[name][label].values())[1][vir_vir[1]]
                vir_1c_occ: float = list(vir_op_vir[name][label].values())[2][vir_vir[0]]
                vir_2c_occ: float = list(vir_op_vir[name][label].values())[2][vir_vir[1]]
                if d_bc:
                    vir_a_occ -= list(occ_op_occ[name][label].values())[0][occ_occ]
                    vir_1b_occ -= list(occ_op_occ[name][label].values())[1][occ_occ]
                    vir_1c_occ -= list(occ_op_occ[name][label].values())[2][occ_occ]
                if d_ad:
                    vir_2b_occ -= list(occ_op_occ[name][label].values())[1][occ_occ]
                    vir_2c_occ -= list(occ_op_occ[name][label].values())[2][occ_occ]
                value[label] = 0.5*lresc_consts[name]*(
                    # <i|A|a>P_{ia,jb}(<b|B|c>-d_{cb}<i|B|j>)P_{ic,jd}<d|C|j>
                (op_a[rotation_i_a[0]]*pp_op_a[0]*vir_1b_occ*pp_op_b[1]*op_c[rotation_j_b[0]+rotations]) +
                    # <i|C|c>P_{ic,jd}(<d|B|a>-d_{ad}<i|B|j>)P_{ia,jb}<b|A|j>
                (op_c[rotation_i_a[1]]*pp_op_b[1]*vir_2b_occ*pp_op_a[0]*op_a[rotation_j_b[1]+rotations]) +
                    # <i|B|a>P_{ia,jb}(<b|A|c>-d_{cb}<i|A|j>)P_{ic,jd}<d|C|j>
                (op_b[rotation_i_a[0]]*pp_op_b[0]*vir_a_occ*pp_op_b[1]*op_c[rotation_j_b[0]+rotations]) +
                    # <i|A|a>P_{ia,jb}(<b|C|c>-d_{cb}<i|C|j>)P_{ic,jd}<d|B|j>
                (op_a[rotation_i_a[0]]*pp_op_a[0]*vir_1c_occ*pp_op_b[1]*op_b[rotation_j_b[0]+rotations]) +
                    # <i|B|c>P_{ic,jd}(<d|C|a>-d_{ad}<i|C|j>)P_{ia,jb}<b|A|j>
                (op_b[rotation_i_a[1]]*pp_op_b[1]*vir_2c_occ*pp_op_a[0]*op_a[rotation_j_b[1]+rotations]) +
                    # <i|C|a>P_{ia,jb}(<b|A|c>-d_{cb}<i|A|j>)P_{ic,jd}<d|B|j>
                (op_c[rotation_i_a[0]]*pp_op_b[0]*vir_a_occ*pp_op_b[1]*op_b[rotation_j_b[0]+rotations]))

                if name in responses.keys():
                    if label in responses[name].keys():
                        value[label] += responses[name][label]
            values[name] = value

    return values


def run_shielding_lresc(wf: wave_function = None, lresc_amounts: list = None,
                                lresc_consts: dict = None, atom: list = None,
                                principal_propagator_approximation: str = None,
                                driver_time: drv_time = None, verbose: int = 0,
                                verbose_integrals: int = -1):
    """
    Calculate all lresc_amounts values
    Args:
        wf (wave_function): Wave function object
        paramagneitc (list): Paramagnetic amount to calculate
        lresc_consts (dict): Dictionary with constants for differente amounts
        atom (list): Atom index to calculate the shielding
        principal_propagator_approximation (str): Approximation
        driver_time (drv_time): Object to manage the time calculation
        verbose (int): Print level
        verbose_integrals (int): Print level for integral module
    """
    start: float = time()

    corrections, responses_amounts, averages = correction_to_calculate(lresc_amounts)
    label_amounts: dict = {"para_nr": paramagnetic_nr,
                    "triplet_lineal_amount": triplet_lineal_responses,
                    "singlet_lineal_amount": singlet_lineal_responses,
                    "triplet_quadratic_amount": triplet_quadratic_responses,
                    "singlet_quadratic_amount": singlet_quadratic_responses,
                    "dia_singlet_lineal": dia_singlet_lineal_responses,
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
    triplet: bool = responses_amounts["triplet_lineal_amount"] or responses_amounts["triplet_quadratic_amount"]
    singlet: bool = (responses_amounts["para_nr"] or responses_amounts["singlet_lineal_amount"]
                    or responses_amounts["dia_singlet_lineal"])
    quadratic_singlet: bool = responses_amounts["singlet_quadratic_amount"]
    if triplet or singlet:
        lineal_pp = drv_principal_propagator(driver_time = driver_time, moe = moe,
                                            n_mo_occ = nocc, n_mo_virt = nvir,
                                            coulomb = coulomb_integrals, exchange = exchange_integrals,
                                            multiplicity_pp = {"singlet": singlet,
                                                                "triplet": triplet},
                                            quadratic = False,
                                            tp_inv = 0, verbose = verbose)
    if quadratic_singlet:
        quadratic_pp = drv_principal_propagator(driver_time = driver_time, moe = moe,
                                            n_mo_occ = nocc, n_mo_virt = nvir,
                                            coulomb = coulomb_integrals, exchange = exchange_integrals,
                                            multiplicity_pp = {"singlet": True,
                                                                "triplet": False},
                                            quadratic = True,
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
                        lresc_consts = lresc_consts, verbose = verbose
                    )
                if activate: atom_av.update(temp_av)
        all_averages[a] = atom_av

        atom_gpvs: dict = {}
        atom_occ_op_occ: dict = {}
        atom_vir_op_vir: dict = {}
        if singlet or triplet or quadratic_singlet:
            for label, activate in responses_amounts.items():
                quadratic: bool = False
                if activate:
                    if "quadratic" in label:
                        quadratic = True
                temp_opo, temp_vpv, temp_gpv = get_gpv(
                        wf = wf, atom = a, driver_time = driver_time, gaugeo = gaugeo,
                        type_correction = label_amounts[label], corrections = corrections,
                        quadratic = quadratic, verbose = verbose, verbose_int = verbose_integrals
                    )
                if activate and quadratic:
                    atom_gpvs.update(temp_gpv)
                    atom_occ_op_occ.update(temp_opo)
                    atom_vir_op_vir.update(temp_vpv)
                elif activate:
                    atom_gpvs.update(temp_gpv)
            gpvs[a] = atom_gpvs
            occ_op_occ[a] = atom_occ_op_occ
            vir_op_vir[a] = atom_vir_op_vir

    all_responses: dict = {}
    lineal_responses: dict = {}
    quadratic_responses: dict = {}
    rotations: int = nocc*nvir
    for i in range(nocc):
        for a in range(nvir):
            s = a + nocc
            for j in range(nocc):
                for b in range(nvir):
                    t = b + nocc

                    if singlet or responses_amounts["triplet_lineal_amount"]:
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
                                    if at not in lineal_responses.keys(): path_responses: dict = {}
                                    else: path_responses: dict = lineal_responses[at]
                                    temp_value: dict = {}
                                    temp_value = get_lineal_response(
                                        responses = path_responses, type_correction = label_amounts[label],
                                        rotations = rotations, corrections = corrections,
                                        lresc_consts = lresc_consts, rotation_i_a = a+i*nvir,
                                        rotation_j_b = b+j*nvir, gpvs = gpvs[at],
                                        pp_left = pp_left, pp_right = pp_right
                                    )
                                if activate: atom_responses.update(temp_value)
                            lineal_responses[at] = atom_responses

                    for c in range(nvir):
                        if not responses_amounts["triplet_quadratic_amount"] and not  responses_amounts["singlet_quadratic_amount"]:
                            break
                        u = c + nocc
                        for d in range(nvir):
                            v = d + nocc

                            for at in atom:
                                atom_responses: dict = {}
                                for label, activate in responses_amounts.items():
                                    if activate and ("triplet_quadratic_amount" == label or "singlet_quadratic_amount" == label):
                                        pp_a: list = [quadratic_pp["singlet"][a+i*nvir,b+j*nvir],
                                                    quadratic_pp["singlet"][c+i*nvir,d+j*nvir]]
                                        if "triplet" in label:
                                            pp_b: list = [lineal_pp["triplet"][a+i*nvir,b+j*nvir],
                                                        lineal_pp["triplet"][c+i*nvir,d+j*nvir]]
                                        else:
                                            pp_b: list = [quadratic_pp["singlet"][a+i*nvir,b+j*nvir],
                                                    quadratic_pp["singlet"][c+i*nvir,d+j*nvir]]
                                        if at not in quadratic_responses.keys(): path_responses: dict = {}
                                        else: path_responses: dict = quadratic_responses[at]
                                        temp_value: dict = {}
                                        temp_value = get_quadratic_response(
                                            responses = path_responses,
                                            type_correction = label_amounts[label],
                                            rotations = rotations, corrections = corrections,
                                            lresc_consts = lresc_consts, rotation_i_a = [a+i*nvir,c+i*nvir],
                                            rotation_j_b = [d+j*nvir, b+j*nvir], occ_occ = j+i*nocc,
                                            vir_vir = [c+b*nvir,a+d*nvir], gpvs = gpvs[at],
                                            occ_op_occ = occ_op_occ[at], vir_op_vir = vir_op_vir[at],
                                            d_bc = b == c, d_ad = a == d, pp_op_a = pp_a, pp_op_b = pp_b
                                        )
                                    if activate: atom_responses.update(temp_value)
                                quadratic_responses[at] = atom_responses

    for at in atom:
        if singlet or responses_amounts["triplet_lineal_amount"]: all_responses[at] = lineal_responses[at]
        if responses_amounts["triplet_quadratic_amount"] and responses_amounts["singlet_quadratic_amount"]:
            if at in all_responses.keys():
                all_responses[at].update(quadratic_responses[at])
            else:
                all_responses[at] = quadratic_responses[at]

    if verbose > 10:
        driver_time.add_name_delta_time(name = "Paramagnetic Calculations", delta_time = (time() - start))

    return all_responses, all_averages

def get_shielding_isotropic(all_responses: dict = None, all_averages: dict = None):

    isotropic_responses: dict = {}
    isotropic_averages: dict = {}

    if all_responses:
        for atom, dict_corrections in  all_responses.items():
            isotropic_atom: dict = None
            for correction, components_values in dict_corrections.items():
                isotropic: dict = {}
                if "fckin" != correction:
                    isotropic[correction] = sum(components_values.values())/3.0
                else:
                    isotropic["fckin"] = sum(components_values.values())
            #
                if not isotropic_atom: isotropic_atom = isotropic
                else: isotropic_atom.update(isotropic)
            #
            isotropic_responses[atom] = isotropic_atom
            if "sdkinxx" in isotropic_responses[atom].keys():
                isotropic_responses[atom]["sdkin"] = (isotropic_responses[atom]["sdkinxx"] +
                                isotropic_responses[atom]["sdkinyy"] +isotropic_responses[atom]["sdkinzz"])
            if "sdbsoxx" in isotropic_responses[atom].keys():
                isotropic_responses[atom]["sdbso"] = (isotropic_responses[atom]["sdbsoxx"] +
                                isotropic_responses[atom]["sdbsoyy"] +isotropic_responses[atom]["sdbsozz"])
            if "lsdsoxx" in isotropic_responses[atom].keys():
                isotropic_responses[atom]["lsdso"] = (isotropic_responses[atom]["lsdsoxx"] +
                                isotropic_responses[atom]["lsdsoyy"] +isotropic_responses[atom]["lsdsozz"])

    if all_averages:
        for atom, dict_corrections in  all_averages.items():
            isotropic_atom: dict = None
            for correction, components_values in dict_corrections.items():
                isotropic: dict = {}
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
                            "lzpso3", "fcsofielzz", "angmomzpsoke3", "psozozke3", "nstcgo3zdw",
                            "nstcgo3zmv", "angmomzfcspinoz", "angmomzpso3massvelo", "angmomzpso3darwin",
                            "nstcgo3z", "psooz3z", "dnske3z"
                        ]
    z_sign: float = -1.0

    if all_responses:
        for atom, dict_corrections in  all_responses.items():
            anisotropic_atom: dict = None
            for correction, components_values in dict_corrections.items():
                anisotropic: dict = {}
                if len(components_values) == 3 and correction not in ["lsdsoxx", "lsdsoyy", "lsdsozz"]:
                    anisotropic[correction] = sum([value if name not in z_component else -value
                                                for name, value in components_values.items()])
                else:
                    anisotropic[correction] = sum(components_values.values())

            #
                if not anisotropic_atom: anisotropic_atom = anisotropic
                else: anisotropic_atom.update(anisotropic)
            #
            anisotropic_responses[atom] = anisotropic_atom
            if "sdkinxx" in anisotropic_responses[atom].keys():
                anisotropic_responses[atom]["sdkin"] = (anisotropic_responses[atom]["sdkinxx"] +
                                anisotropic_responses[atom]["sdkinyy"] + z_sign*anisotropic_responses[atom]["sdkinzz"])
            if "sdbsoxx" in anisotropic_responses[atom].keys():
                anisotropic_responses[atom]["sdbso"] = (anisotropic_responses[atom]["sdbsoxx"] +
                                anisotropic_responses[atom]["sdbsoyy"] + z_sign*anisotropic_responses[atom]["sdbsozz"])
            if "lsdsoxx" in anisotropic_responses[atom].keys():
                anisotropic_responses[atom]["lsdso"] = (anisotropic_responses[atom]["lsdsoxx"] +
                                anisotropic_responses[atom]["lsdsoyy"] + z_sign*anisotropic_responses[atom]["lsdsozz"])

    if all_averages:
        for atom, dict_corrections in  all_averages.items():
            anisotropic_atom: dict = None
            for correction, components_values in dict_corrections.items():
                anisotropic: dict = {}
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
    if "fckin" in responses.keys():
        fckin: float = responses["fckin"]
        ani_fckin: float = ani_responses["fckin"]
    else:
        fckin: str = "No Calc."
        ani_fckin: str = "No Calc."
    if "sdkin" in responses.keys():
        sdkin: float = responses["sdkin"]
        ani_sdkin: float = ani_responses["sdkin"]
    else:
        sdkin: str = "No Calc."
        ani_sdkin: str = "No Calc."
    if "fcbso" in responses.keys():
        fcbso: float = responses["fcbso"]
        ani_fcbso: float = ani_responses["fcbso"]
    else:
        fcbso: str = "No Calc."
        ani_fcbso: str = "No Calc."
    if "sdbso" in responses.keys():
        sdbso: float = responses["sdbso"]
        ani_sdbso: float = ani_responses["sdbso"]
    else:
        sdbso: str = "No Calc."
        ani_sdbso: str = "No Calc."
    print(("<<FC; K>>".center(20) + "   " + "<<SD; K>>".center(20) + "    " +
            "<<FC; BSO>>".center(20) + "   " + "<<SD; BSO>>".center(20)).center(89))
    print(("|" + "¯"*18 + "|" + "   " + "|" + "¯"*18 + "|" + "   " +
            "|" + "¯"*18 + "|" + "   " + "|" + "¯"*18 + "|").center(89))
    print(("|" + f"{fckin:.6f}".center(18) + "|" + "   " + "|" + f"{sdkin:.6f}".center(18) + "|"
            + "   " + "|" + f"{fcbso:.6f}".center(18) + "|" + "   " + "|" + f"{sdbso:.6f}".center(18) + "|").center(89))
    print(("|" + f"{ani_fckin:.6f}".center(18) + "|" + "   " + "|" + f"{ani_sdkin:.6f}".center(18) + "|"
            + "   " + "|" + f"{ani_fcbso:.6f}".center(18) + "|" + "   " + "|" + f"{ani_sdbso:.6f}".center(18) + "|").center(89))
    print(("¯"*18 + "     " + "¯"*18 + "     " + "¯"*18 + "     " + "¯"*18).center(89))

def print_sinlget_quadratic(responses: dict = None, ani_responses: dict = None):
    if "lpsomv" in responses.keys():
        lpsomv: float = responses["lpsomv"]
        ani_lpsomv: float = ani_responses["lpsomv"]
    else:
        lpsomv: str = "No Calc."
        ani_lpsomv: str = "No Calc."
    if "lpsodw" in responses.keys():
        lpsodw: float = responses["lpsodw"]
        ani_lpsodw: float = ani_responses["lpsodw"]
    else:
        lpsodw: str = "No Calc."
        ani_lpsodw: str = "No Calc."

    print(("<<L; PSO, Mv>>".center(20) + "   " + "<<L; PSO, Dw>>".center(20)).center(89))
    print(("|" + "¯"*18 + "|" + "   " + "|" + "¯"*18 + "|").center(89))
    print(("|" + f"{lpsomv:.6f}".center(18) + "|" + "   " + "|" + f"{lpsodw:.6f}".center(18) + "|").center(89))
    print(("|" + f"{ani_lpsomv:.6f}".center(18) + "|" + "   " + "|" + f"{ani_lpsodw:.6f}".center(18) + "|").center(89))
    print(("¯"*18 + "     " + "¯"*18).center(89))

def print_triplet_quadratic(responses: dict = None, ani_responses: dict = None):
    if "lfcso" in responses.keys():
        lfcso: float = responses["lfcso"]
        ani_lfcso: float = ani_responses["lfcso"]
    else:
        lfcso: str = "No Calc."
        ani_lfcso: str = "No Calc."
    if "lsdso" in responses.keys():
        lsdso: float = responses["lsdso"]
        ani_lsdso: float = ani_responses["lsdso"]
    else:
        lsdso: str = "No Calc."
        ani_lsdso: str = "No Calc."

    print(("<<L; FC, SO>>".center(20) + "   " + "<<L; SD, SO>>".center(20)).center(89))
    print(("|" + "¯"*18 + "|" + "   " + "|" + "¯"*18 + "|").center(89))
    print(("|" + f"{lfcso:.6f}".center(18) + "|" + "   " + "|" + f"{lsdso:.6f}".center(18) + "|").center(89))
    print(("|" + f"{ani_lfcso:.6f}".center(18) + "|" + "   "
            + "|" + f"{ani_lsdso:.6f}".center(18) + "|").center(89))
    print(("¯"*18 + "     " + "¯"*18).center(89))

def print_dia_sinlget_lineal(responses: dict = None, ani_responses: dict = None):
    if "a2mv" in responses.keys():
        a2mv: float = responses["a2mv"]
        ani_a2mv: float = ani_responses["a2mv"]
    else:
        a2mv: str = "No Calc."
        ani_a2mv: str = "No Calc."
    if "a2dw" in responses.keys():
        a2dw: float = responses["a2dw"]
        ani_a2dw: float = ani_responses["a2dw"]
    else:
        a2dw: str = "No Calc."
        ani_a2dw: str = "No Calc."

    print(("<<A²; Mv>>".center(20) + "   " + "<<A²; Dw>>".center(20)).center(89))
    print(("|" + "¯"*18 + "|" + "   " + "|" + "¯"*18 + "|").center(89))
    print(("|" + f"{a2mv:.6f}".center(18) + "|" + "   " + "|" + f"{a2dw:.6f}".center(18) + "|").center(89))
    print(("|" + f"{ani_a2mv:.6f}".center(18) + "|" + "   "
            + "|" + f"{ani_a2dw:.6f}".center(18) + "|").center(89))
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


    print(("<FC>".center(20) + "   " + "<PSOOZ>".center(20) + "   " + "<DNSKE>".center(20)).center(89))
    print(("|" + "¯"*18 + "|" + "   " + "|" + "¯"*18 + "|" + "   " + "|" + "¯"*18 + "|").center(89))
    print(("|" + f"{fc:.6f}".center(18) + "|" + "   " + "|" + f"{psooz:.6f}".center(18) + "|"
            + "   " + "|" + f"{dnske:.6f}".center(18) + "|").center(89))
    print(("|" + f"{ani_fc:.6f}".center(18) + "|" + "   " + "|" + f"{ani_psooz:.6f}".center(18) + "|"
            + "   " + "|" + f"{ani_dnske:.6f}".center(18) + "|").center(89))
    print(("¯"*18 + "     " + "¯"*18  + "     " + "¯"*18).center(89))

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

    print(f"Isotropic Values [ppm]".center(89), "\n", "and".center(89), "Anisotropic Values [ppm]".center(89))
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
        print_sinlget_quadratic(isotropic_responses[atom],anisotropic_responses[atom])
        print_triplet_quadratic(isotropic_responses[atom],anisotropic_responses[atom])

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
        print_dia_sinlget_lineal(isotropic_responses[atom],anisotropic_responses[atom])
        print_dia_averages(isotropic_averages[atom], anisotropic_averages[atom])

    # Brief results
        # -- isotropic
        lresc_para = sum([correction for label, correction in isotropic_responses[atom].items()
                            if label not in ["paranr", "a2mv", "a2dw"]])
        lresc_dia = (sum([correction for label, correction in isotropic_responses[atom].items()
                            if label in ["a2mv", "a2dw"]])
                        + sum([correction for label, correction in isotropic_averages[atom].items()
                        if label not in ["dianr"]]))
        ligand_correction = sum([correction for label, correction in isotropic_responses[atom].items()
                            if label in ["lfcso", "lsdso", "lpsodw", "lpsomv", "lpsokin", "lkinpso"]])
        core_correction = (sum([correction for label, correction in isotropic_averages[atom].items()
                            if label not in ["dianr"]])
                            + sum([correction for label, correction in isotropic_responses[atom].items()
                        if label not in ["lfcso", "lsdso", "lpsodw", "lpsomv", "lpsokin", "lkinpso"]]))
        # -- anisotropic
        ani_lresc_para = sum([correction for label, correction in anisotropic_responses[atom].items()
                            if label not in ["paranr", "a2mv", "a2dw"]])
        ani_lresc_dia = (sum([correction for label, correction in anisotropic_responses[atom].items()
                            if label in ["a2mv", "a2dw"]])
                        + sum([correction for label, correction in anisotropic_averages[atom].items()
                        if label not in ["dianr"]]))
        ani_ligand_correction = sum([correction for label, correction in anisotropic_responses[atom].items()
                            if label in ["lfcso", "lsdso", "lpsodw", "lpsomv", "lpsokin", "lkinpso"]])
        ani_core_correction = (sum([correction for label, correction in anisotropic_averages[atom].items()
                            if label not in ["dianr"]])
                            + sum([correction for label, correction in anisotropic_responses[atom].items()
                        if label not in ["lfcso", "lsdso", "lpsodw", "lpsomv", "lpsokin", "lkinpso"]]))

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
