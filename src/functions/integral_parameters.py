from libf import *

def integral_1b_parameters(atoms_number: int = None, integral_name: str = None, verbose: int = 0):
    """
    Build dictionary with information neccesary to calculate of the
    integrals

    Args:
        wf (wave_function): Wave Function object
        integral_name (str): Integral name.
        verbose (int): Print level.
    """
    # Default values for integrals properties ######################################################
    r_gauge: list = [0.0,0.0,0.0]
    r_dipole: list = [0.0,0.0,0.0]
    atoms_number: int = atoms_number
    atoms: list = [atom for atom in range(atoms_number)]
    if integral_name.lower() == "laplacian":
        magnetic_components: list = [0,1,2,3,4,5]
    else:
        magnetic_components: list = [0,1,2]
    spatial_symmetries: list = [symmetry for symmetry in range(atoms_number*3)]
    ###############################################################################################
    property_split = integral_name.lower().split()

    if len(property_split) > 1:

        # label symmetry magnetic
        if spatial_symmetry[property_split[0]] == 1 and magnetic[property_split[0]] == 1:
            spatial_symmetries = [int(property_split[1]) - 1]
            if property_split[2].isnumeric():
                magnetic_components = [int(property_split[2]) - 1]
            else:
                if property_split.lower() == "laplacian":
                    magnetic_components = [i for i, x in spatial_components.items() if x == property_split[2]]
                else:
                    magnetic_components = [i for i, b in magnetic_axes.items() if b == property_split[2]]

        # label symmetry
        elif spatial_symmetry[property_split[0]] == 1  and magnetic[property_split[0]] == 0:
            spatial_symmetries = [int(property_split[1]) - 1]

        # label magnetic
        elif magnetic[property_split[0]] == 1 and spatial_symmetry[property_split[0]] == 0:
            if property_split[1].isnumeric():
                magnetic_components = [int(property_split[1]) - 1]
            else:
                if property_split.lower() == "laplacian":
                    magnetic_components = [i for i, x in spatial_components.items() if x == property_split[2]]
                else:
                    magnetic_components = [i for i, b in magnetic_axes.items() if b == property_split[1]]

        elif property_split[0] in ["fc", "nucpot"]:
            atoms = [int(property_split[1]) - 1]

    return r_gauge, r_dipole, magnetic_components, spatial_symmetries, atoms
    #return all_response, parameters