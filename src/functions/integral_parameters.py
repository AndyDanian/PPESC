from libf import *

def integral_1b_parameters(wf: wave_function = None, property: str = None, verbose: int = 0):
    """
    Build dictionary with information neccesary to calculate of the
    integrals

    Args:
        wf (wave_function): Wave Function object
        property (str): Property label.
        verbose (int): Print level.
    """
    property_split = property.lower().split()
    if property.lower() == "fc" or "fc" in property_split:
        if property.lower() == "fc":
            parameters = None #Calculate Fermi--contact on all atoms
            all_response = True
        else:
            parameters = {"atoms": [int(property_split[1])]}
            all_response = False

    return all_response, parameters