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
            parameters = {"atoms": [a for a in range(wf.atom_number)]}
        else:
            parameters = {"atoms": [int(property_split)]}

    return parameters