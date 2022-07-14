from libl import *

def correction_to_calculate(paramagnetic: list = None):
    """
    Turn on or off paramagnetic corrections
    Args:
        paramagneitc (list): Paramagnetic amount to calculate
    """
    triplet_lineal_amount: bool = False
    singlet_lineal_amount: bool = False
    triplet_quadratic_amount: bool = False
    singlet_quadratic_amount: bool = False
    if not paramagnetic:
        return {
        # Lineal -- Triplet
        "fckin": True, "sdkinxx": True, "sdkinyy": True, "sdkinzz": True,
        "fcbso": True, "sdbsoxx": True, "sdbsoyy": True, "sdbsozz": True,
        # Lineal -- Singlet
        "lpsokin": True, "lkinpso": True,
        # Quadratic -- Triplet
        "lfcso" : True,
        "lsdsoxx": True, "lsdsoyy": True, "lsdsozz": True,
        # Quadratic -- Singlet
        "lpsodw": True, "lpsomv": True
        }, True, True, True, True
    else:
        if ("fckin" in paramagnetic or "sdkin" in paramagnetic or
            "fcbso" in paramagnetic or "sdbso" in paramagnetic): triplet_lineal_amount = True
        if "lpsokin" in paramagnetic or "lkinpso" in paramagnetic: singlet_lineal_amount = True
        if ("lfcso" in paramagnetic or "lsdsoxx" in paramagnetic or
            "lsdsoyy" in paramagnetic or "lsdsozz" in paramagnetic): triplet_quadratic_amount = True
        if ("lpsodw" in paramagnetic or "lpsomv" in paramagnetic): singlet_quadratic_amount = True
        return {
        # Lineal -- Triplet
        "fckin": "fckin" in paramagnetic, "sdkin": "sdkin" in paramagnetic,
        "fcbso": "fcbso" in paramagnetic, "sdbso": "sdbso" in paramagnetic,
        # Lineal -- Singlet
        "lpsokin": "lpsokin" in paramagnetic, "lkinpso": "lkinpso" in paramagnetic,
        # Quadratic -- Triplet
        "lfcso" : "lfcso" in paramagnetic,
        "lsdsoxx": "lsdsoxx" in paramagnetic, "lsdsoyy": "lsdsoyy" in paramagnetic, "lsdsozz": "lsdsozz" in paramagnetic,
        # Quadratic -- Singlet
        "lpsodw": "lpsodw" in paramagnetic, "lpsomv": "lpsomv" in paramagnetic
        }, triplet_lineal_amount, singlet_lineal_amount, triplet_quadratic_amount, singlet_quadratic_amount

def get_responses(rps: response = None, atom: int = None, corrections: dict = None,
                    lresc_consts: dict = None, type_correction: dict = None,
                    principal_propagator_approximation: str = None,
                    gaugeo: list = None, verbose_responses: int = -1):
    """_summary_

    Args:
        rps (response): Response object
        atom (int): Atom index
        corrections (dict): Dictionary with corrections to calculate
        lresc_consts (dict): Dictionary with constants for differente amounts
        type_correction (dict): Type of correction
        principal_propagator_approximation (str): Approximation for principal propagator
        gaugeo (list): Gauge
        verbose_responses (int): Print level for response
    """
    responses_corrections: dict = {}
    for name, amount in type_correction.items():
        if corrections[name]:
            response_correction: dict = {}
            for label, response in amount.items():
                response_correction[label] = (lresc_consts[name] *
                                        list(rps.drv_reponse_calculation(principal_propagator_approximation = principal_propagator_approximation,
                                                properties = [response(atom)], gaugeo = gaugeo, verbose = verbose_responses).values())[0]
                                        )
            print(" corrections  ", response_correction)
            responses_corrections[name] = response_correction
    return responses_corrections

# def print_paramagnetic():
#     """
#     Print the paramagnetic amounts
#     """

def run_shielding_paramagnetic(wf: wave_function = None, paramagnetic: list = None,
                                lresc_consts: dict = None, atom: list = None,
                                principal_propagator_approximation: str = None,
                                driver_time: drv_time = None, verbose: int = 0,
                                verbose_responses: int = -1):
    """
    Calculate all paramagnetic values
    Args:
        wf (wave_function): Wave function object
        paramagneitc (list): Paramagnetic amount to calculate
        lresc_consts (dict): Dictionary with constants for differente amounts
        atom (list): Atom index to calculate the shielding
        principal_propagator_approximation (str): Approximation
        driver_time (drv_time): Object to manage the time calculation
        verbose (int): Print level
        verbose_responses (int): Print level for responses module
    """
    start: float = time()

    corrections, triplet_lineal, singlet_lineal, triplet_quadratic, singlet_quadratic = correction_to_calculate(paramagnetic)

    for a in atom:
        rps = response(wf = wf)
        gaugeo: list = wf.coordinates[a]
        # -- Lineal Response
        # - Triplet
        if triplet_lineal:
            triplet_lineal_corrections = get_responses(rps = rps, corrections = corrections, atom = a, lresc_consts = lresc_consts,
                                            type_correction = triplet_lineal_responses,
                                            principal_propagator_approximation = principal_propagator_approximation,
                                            gaugeo = gaugeo, verbose_responses = verbose_responses)
        # FcKin
        # if corrections["fckin"]:
        #     responses_fckin: dict = {}
        #     responses_fckin["fckin"] = (lresc_consts["fckin"] *
        #                                 list(rps.drv_reponse_calculation(principal_propagator_approximation = principal_propagator_approximation,
        #                                         properties = [fckin["fckin"](a)], gaugeo = gaugeo, verbose = verbose_responses).values())[0]
        #                                 )
        #     print("fckin ",responses_fckin)
        # # SdKin
        # if corrections["sdkinxx"]:
        #     responses_sdkinxx: dict = {}
        #     for name, value in sdkinxx.items():
        #         responses_sdkinxx[name] = (lresc_consts["sdkin"] *
        #                                     list(rps.drv_reponse_calculation(principal_propagator_approximation = principal_propagator_approximation,
        #                                         properties = [value(a)], gaugeo = gaugeo, verbose = verbose_responses).values())[0])
        # if corrections["sdkinyy"]:
        #     responses_sdkinyy: dict = {}
        #     for name, value in sdkinyy.items():
        #         responses_sdkinyy[name] = (lresc_consts["sdkin"] *
        #                                     list(rps.drv_reponse_calculation(principal_propagator_approximation = principal_propagator_approximation,
        #                                         properties = [value(a)], gaugeo = gaugeo, verbose = verbose_responses).values())[0])
        # if corrections["sdkinzz"]:
        #     responses_sdkinzz: dict = {}
        #     for name, value in sdkinzz.items():
        #         responses_sdkinzz[name] = (lresc_consts["sdkin"] *
        #                                     list(rps.drv_reponse_calculation(principal_propagator_approximation = principal_propagator_approximation,
        #                                         properties = [value(a)], gaugeo = gaugeo, verbose = verbose_responses).values())[0])
        #     print("sdkin ",(sum(list(responses_sdkinxx.values())) +
        #                     sum(list(responses_sdkinyy.values())) +
        #                     sum(list(responses_sdkinzz.values())))/3.0)
        # # FcBso
        # if corrections["fcbso"]:
        #     responses_fcbso: dict = {}
        #     for name, value in fcbso.items():
        #         responses_fcbso[name] = (lresc_consts["fcbso"] *
        #                                 list(rps.drv_reponse_calculation(principal_propagator_approximation = principal_propagator_approximation,
        #                                         properties = [value(a)], gaugeo = gaugeo, verbose = verbose_responses).values())[0])
        #     print("fcbso ",sum(list(responses_fcbso.values()))/3.0)
        # # SdBso
        # if corrections["sdbsoxx"]:
        #     responses_sdbsoxx: dict = {}
        #     for name, value in sdbsoxx.items():
        #         responses_sdbsoxx[name] = (lresc_consts["sdbso"] *
        #                                     list(rps.drv_reponse_calculation(principal_propagator_approximation = principal_propagator_approximation,
        #                                         properties = [value(a)], gaugeo = gaugeo, verbose = verbose_responses).values())[0])
        # if corrections["sdbsoyy"]:
        #     responses_sdbsoyy: dict = {}
        #     for name, value in sdbsoyy.items():
        #         responses_sdbsoyy[name] = (lresc_consts["sdbso"] *
        #                                     list(rps.drv_reponse_calculation(principal_propagator_approximation = principal_propagator_approximation,
        #                                         properties = [value(a)], gaugeo = gaugeo, verbose = verbose_responses).values())[0])
        # if corrections["sdbsozz"]:
        #     responses_sdbsozz: dict = {}
        #     for name, value in sdbsozz.items():
        #         responses_sdbsozz[name] = (lresc_consts["sdbso"] *
        #                                     list(rps.drv_reponse_calculation(principal_propagator_approximation = principal_propagator_approximation,
        #                                         properties = [value(a)], gaugeo = gaugeo, verbose = verbose_responses).values())[0])
        #     print("sdbso ",(sum(list(responses_sdbsoxx.values())) +
        #                     sum(list(responses_sdbsoyy.values())) +
        #                     sum(list(responses_sdbsozz.values())))/3.0)

        # - Singlet
        if singlet_lineal:
            singlet_lineal_corrections = get_responses(rps = rps, corrections = corrections, atom = a, lresc_consts = lresc_consts,
                                                type_correction = singlet_lineal_responses,
                                                principal_propagator_approximation = principal_propagator_approximation,
                                                gaugeo = gaugeo, verbose_responses = verbose_responses)
        # # L-PsoKin
        # if corrections["lpsokin"]:
        #     responses_lpsokin: dict = {}
        #     for name, value in lpsokin.items():
        #         responses_lpsokin[name] = (lresc_consts["lpsokin"] *
        #                                     list(rps.drv_reponse_calculation(principal_propagator_approximation = principal_propagator_approximation,
        #                                         properties = [value(a)], gaugeo = gaugeo, verbose = verbose_responses).values())[0])
        #     print("lpsokin ", sum(list(responses_lpsokin.values()))/3.0)
        # if corrections["lkinpso"]:
        #     responses_lkinpso: dict = {}
        #     for name, value in lkinpso.items():
        #         responses_lkinpso[name] = (lresc_consts["lkinpso"] *
        #                                     list(rps.drv_reponse_calculation(principal_propagator_approximation = principal_propagator_approximation,
        #                                         properties = [value(a)], gaugeo = gaugeo, verbose = verbose_responses).values())[0])
        #     print("lkinpso ", sum(list(responses_lkinpso.values()))/3.0)
        # -- Quadratic Responses
        # - Triplet
        if triplet_quadratic:
            triplet_quadratic_corrections = get_responses(rps = rps, corrections = corrections, atom = a, lresc_consts = lresc_consts,
                                                type_correction = triplet_quadratic_responses,
                                                principal_propagator_approximation = principal_propagator_approximation,
                                                gaugeo = gaugeo, verbose_responses = verbose_responses)
        # # LFcSO
        # if corrections["lfcso"]:
        #     responses_lfcso: dict = {}
        #     for name, value in lfcso.items():
        #         responses_lfcso[name] = (lresc_consts["lfcso"] *
        #                                 list(rps.drv_reponse_calculation(principal_propagator_approximation = principal_propagator_approximation,
        #                                         properties = [value(a)], gaugeo = gaugeo, verbose = verbose_responses).values())[0])
        #     print("lfcso ", sum(list(responses_lfcso.values()))/3.0)
        # #LSdSOxx
        # if corrections["lsdsoxx"]:
        #     responses_lsdsoxx: dict = {}
        #     for name, value in lsdsoxx.items():
        #         responses_lsdsoxx[name] = (lresc_consts["lsdso"] *
        #                                     list(rps.drv_reponse_calculation(principal_propagator_approximation = principal_propagator_approximation,
        #                                         properties = [value(a)], gaugeo = gaugeo, verbose = verbose_responses).values())[0])
        # #LSdSOxx
        # if corrections["lsdsoyy"]:
        #     responses_lsdsoyy: dict = {}
        #     for name, value in lsdsoyy.items():
        #         responses_lsdsoyy[name] = (lresc_consts["lsdso"] *
        #                                     list(rps.drv_reponse_calculation(principal_propagator_approximation = principal_propagator_approximation,
        #                                         properties = [value(a)], gaugeo = gaugeo, verbose = verbose_responses).values())[0])
        # #LSdSOxx
        # if corrections["lsdsozz"]:
        #     responses_lsdsozz: dict = {}
        #     for name, value in lsdsozz.items():
        #         responses_lsdsozz[name] = (lresc_consts["lsdso"] *
        #                                     list(rps.drv_reponse_calculation(principal_propagator_approximation = principal_propagator_approximation,
        #                                         properties = [value(a)], gaugeo = gaugeo, verbose = verbose_responses).values())[0])
        #     print("lsdso ",(sum(list(responses_lsdsoxx.values())) +
        #                     sum(list(responses_lsdsoyy.values())) +
        #                     sum(list(responses_lsdsozz.values())))/3.0)
        # - Singlet
        if singlet_quadratic:
            singlet_quadratic_corrections = get_responses(rps = rps, corrections = corrections, atom = a, lresc_consts = lresc_consts,
                                                type_correction = singlet_quadratic_responses,
                                                principal_propagator_approximation = principal_propagator_approximation,
                                                gaugeo = gaugeo, verbose_responses = verbose_responses)
        # # LpsoMv
        # if corrections["lpsomv"]:
        #     responses_lpsomv: dict = {}
        #     for name, value in lpsomv.items():
        #         responses_lpsomv[name] = (lresc_consts["lpsomv"] *
        #                                 list(rps.drv_reponse_calculation(principal_propagator_approximation = principal_propagator_approximation,
        #                                         properties = [value(a)], gaugeo = gaugeo, verbose = verbose_responses).values())[0])
        #     print("lpsomv ", sum(list(responses_lpsomv.values()))/3.0)
        # # LpsoDw
        # if corrections["lpsodw"]:
        #     responses_lpsodw: dict = {}
        #     for name, value in lpsodw.items():
        #         responses_lpsodw[name] = (lresc_consts["lpsodw"] *
        #                                 list(rps.drv_reponse_calculation(principal_propagator_approximation = principal_propagator_approximation,
        #                                         properties = [value(a)], gaugeo = gaugeo, verbose = verbose_responses).values())[0])
        #     print("lpsodw ", sum(list(responses_lpsodw.values()))/3.0)

        #print_paramagnetic(triplet_lineal_corrections, singlet_lineal_corrections, triplet_quadratic_corrections, singlet_quadratic_corrections)

    if verbose > 10:
        driver_time.add_name_delta_time(name = "Paramagnetic Calculations", delta_time = (time() - start))