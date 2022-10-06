from libf import *


def integral_1b_parameters(
    atoms_number: int,
    integral_name: str,
    gauge: list[float],
    dipole: list[float],
    verbose: int = 0,
):
    """
    Build dictionary with information neccesary to calculate of the
    integrals

    Args:
        wf (wave_function): Wave Function object
        integral_name (str): Integral name.
        gauge (list): Gauge origen
        dipole (list): Dipole origen
        verbose (int): Print level.
    """
    # Default values for integrals properties ######################################################
    if gauge:
        r_gauge: list[float] = gauge
    else:
        r_gauge = [0.0, 0.0, 0.0]

    if dipole:
        r_dipole: list[float] = dipole
    else:
        r_dipole = [0.0, 0.0, 0.0]

    atoms: list = [atom for atom in range(atoms_number)]

    if integral_name.lower() == "laplacian":
        magnetic_components: list[int] = [0, 1, 2, 3, 4, 5]
    else:
        magnetic_components = [0, 1, 2]

    spatial_symmetries: list[int] = [symmetry for symmetry in range(atoms_number * 3)]
    ###############################################################################################
    property_split: list[str] = integral_name.lower().split()

    if len(property_split) > 1:
        # label symmetry magnetic
        if (
            spatial_symmetry[property_split[0]] == 1
            and magnetic[property_split[0]] == 1
        ):
            spatial_symmetries = [int(property_split[1]) - 1]
            if property_split[2].isnumeric():
                magnetic_components = [int(property_split[2]) - 1]
            else:
                magnetic_components = [
                    i for i, b in magnetic_axes.items() if b == property_split[2]
                ]
        # label symmetry
        elif (
            spatial_symmetry[property_split[0]] == 1
            and magnetic[property_split[0]] == 0
        ):
            spatial_symmetries = [int(property_split[1]) - 1]

        # label magnetic
        elif (
            magnetic[property_split[0]] == 1
            and spatial_symmetry[property_split[0]] == 0
        ):
            if property_split[1].isnumeric():
                magnetic_components = [int(property_split[1]) - 1]
            else:
                if "laplacian" in property_split:
                    magnetic_components = [
                        property_split[1] if property_split[1].isnumeric() else i
                        for i, x in spatial_components.items()
                        if x == property_split[1]
                    ]
                else:
                    magnetic_components = [
                        i for i, b in magnetic_axes.items() if b == property_split[1]
                    ]

        elif property_split[0] in ["fc", "nucpot"]:
            atoms = [int(property_split[1]) - 1]

    return r_gauge, r_dipole, magnetic_components, spatial_symmetries, atoms
