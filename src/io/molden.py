from os import read
from libsrc import *
import numpy as np

def convert_to_float(variable_name: None, value: None):
    """
    Convert numeric string o int to float or return None
    when is a special character like ****
    """
    try:
        return float(value)
    except:
        print(f"**** Warning\n\n\
        Value {value} assigned to {variable_name} isn't a float.\n\
        This can produce errors in the calcule when is neccesary {variable_name}.\n")
        return None

def validate_file(file_molden):
    """
    File validation

    Args:
        file_molden (file): Molden's File
    """
    gto_section = mo_section = atoms_section = primitive_type = False

    with open(file_molden, "r") as f:
        content = f.readlines()
        for element in content:
            if "[GTO]" in element:
                gto_section = True
            elif "[MO]" in element:
                mo_section = True
            elif "[Atoms]" in element:
                atoms_section = True
            elif "[5D" in element or "[6D" in element:
                primitive_type = True

            if gto_section and mo_section and atoms_section:
                return True
                break

        if not gto_section:
            ValueError("*** Error\n There are not bases set informations.")
        elif not mo_section:
            ValueError("*** Error\n There are not MO Coefficientes.")
        elif not atoms_section:
            ValueError("*** Error\n There are not Atomic information.")
        elif not primitive_type:
            ValueError(
                "*** Error\n\
            There is not information about primitive type\
            \n spherical or cartessian ([6D..] or [5D..])."
            )


def read_molden(file_molden: str = None, drv_scratch: scratch = None, verbose=21):
    """[summary]
    Read .molden generate with DALTON the most general possible
    Args:
        file_molden (string): Name .molden
    Return:
        Z_l_xyz (array)         : Atomic number, label and coordinate
        Primitive_Atoms (array) : Pre-exponential and contraction coefficients
        Amount_TP_ML (array)    : Among by each type of primitive of each atom
        TP (array)              : Symbol of type of primitive
        Occ_MO   (array)        : Occupation of each of the OMs
        Coeff_MO (array)        : Molecular coefficcients by each MO
    """
    label_title: str = "Reading Wave Function From MOLDEN File"
    drv_scratch.write_output(label_title.center(80))
    drv_scratch.write_output(("-"*len(label_title)).center(80))

    # Verification of the file
    if validate_file(file_molden):
        drv_scratch.write_output("Molden's file is valid")

    # Number of lines where are the information
    content = open(file_molden, "r")
    datafile = content.readlines()
    amount_strings = 0

    # label, index, charge, and coordinates
    numbers_atoms = 0
    A_Bohr = 1.0
    stop_atoms = False
    molecule_index = 0
    molecule_index_before = 0
    molecule_information = []
    l_i_q_xyz = []

    # Primitive type and amount, and exponents
    # NOTE: Only bases set uncontracteds
    stop_gto = False
    t_a_exp = []
    line_gto = 2
    stop_verification_contraction = False
    total_primitives = 0
    # -- Molden Format https://www3.cmbi.umcn.nl/molden/molden_format.html
    l = {
        "cartessian": {"s": 1, "p": 3, "d": 6, "f": 10, "g": 15, "h": 21, "i": 28},
        "spherical": {"s": 1, "p": 3, "d": 5, "f": 7, "g": 9, "h": 11, "i": 13},
    }
    spatial_primitive = "cartessian"
    #
    spatial_primitive_find = False
    count_spatial = 0

    # MO Coefficients
    stop_mo = False
    count_mo = 1
    mo = []
    num_mo = 1
    end_mo = ["S", "[", "E"]

    for number_line, element in enumerate(datafile):
        if "[TITLE]" in element:
            drv_scratch.write_output("Title into Molden: ")
            drv_scratch.write_output(datafile[number_line + 1])

        while "[Atoms]" in element and stop_atoms == False:
            # Read atomic type and coordinates
            #
            # [Atoms] Au/Angstrom
            # Atom Symbol  Atomic Number  Charge  x  y  z
            # Atom Symbol  Atomic Number  Charge  x  y  z
            #
            if datafile[number_line].split()[1] != "AU":
                A_Bohr = 1.0 / 0.52917721

            # Geometry last line
            # it is neccesary because there are sometimes where the next string
            # have six string
            no_float = False
            if len(datafile[number_line + 1 + numbers_atoms].split()) >= 3:
                third_string = datafile[number_line + 1 + numbers_atoms].split()[3]
                no_float = third_string.isalpha()

            if (
                # Amount the strings, molden format is 6 to geometry
                (amount_strings != len(datafile[number_line + 1 + numbers_atoms].split())
                or
                # If the third string in the line is not float then stop
                no_float)
                and numbers_atoms > 0
            ):
                stop_atoms = True
                if len(l_i_q_xyz) == numbers_atoms:
                    print(
                        "*** Warning VER ESTO MEJORAR\n\
                    Each atom is a molecule. You can group atoms put on the same index\n\
                    into molden file in section [Atoms], second column. Example: \n\n\
                    O 1 8 0.0 0.0 0.0\n\
                    H 1 1 0.0 0.0 0.0\n\
                    H 1 1 0.0 0.0 0.0\n"
                    )
            else:
                molecule_information = [
                    datafile[number_line + 1 + numbers_atoms].split()[0]
                    + " "
                    + datafile[number_line + 1 + numbers_atoms].split()[2]
                    + " "
                    + str(
                        convert_to_float(
                            "x coordinate", datafile[number_line + 1 + numbers_atoms].split()[
                                3
                            ]
                        )
                        * A_Bohr
                    )
                    + " "
                    + str(
                        convert_to_float(
                            "y coordinate", datafile[number_line + 1 + numbers_atoms].split()[
                                4
                            ]
                        )
                        * A_Bohr
                    )
                    + " "
                    + str(
                        convert_to_float(
                            "z coordinate",
                            datafile[number_line + 1 + numbers_atoms].split()[
                                5
                            ]
                        )
                        * A_Bohr
                    )
                ]

                molecule_index = int(
                    datafile[number_line + 1 + numbers_atoms].split()[1]
                )
                if (
                    numbers_atoms > 0
                    and molecule_index != molecule_index_before
                ):
                    l_i_q_xyz.append(molecule_information)

                elif numbers_atoms == 0:
                    l_i_q_xyz.append(molecule_information)

                elif (
                    numbers_atoms > 0
                    and molecule_index == molecule_index_before
                ):
                    l_i_q_xyz[numbers_atoms - 1] += molecule_information

                molecule_index_before = molecule_index

                amount_strings = len(
                    datafile[number_line + 1 + numbers_atoms].split()
                )
                numbers_atoms += 1

        # Cartessian or Spherical
        while "5D" in element and spatial_primitive_find == False:
            pass
            spatial_primitive_find = True
        else:
            if spatial_primitive_find == False:
                spatial_primitive = "spherical"
                spatial_primitive_find = True

        while "[GTO]" in element and stop_gto == False:

            while spatial_primitive_find == False:
                datafile[number_line + count_spatial]
                if "5D" in datafile[number_line + count_spatial]:
                    pass
                    spatial_primitive_find = True
                else:
                    spatial_primitive = "spherical"
                    spatial_primitive_find = True
                count_spatial += 1

            # GTO information
            #
            # [GTO]
            # atomic index  0
            # primitive type   primitive amount   1.00
            # primitive exponents   contraction coefficients
            #
            primitive = {}
            stop_read_exponents = False
            while stop_read_exponents == False:
                # Primitive type
                primitive_type = datafile[number_line + line_gto].split()[0]
                # Primitive amount
                primitive_amount = int(
                    datafile[number_line + line_gto].split()[1]
                )
                line_gto += 1
                # Exponents
                exponents = [
                    float(datafile[number_line + line_gto + i].split()[0])
                    for i in range(primitive_amount)
                ]

                # Contraction coefficients
                contraction_coefficients = [
                    float(datafile[number_line + line_gto + i].split()[1])
                    for i in range(primitive_amount)
                ]

                if not stop_verification_contraction:
                    for coefficient in contraction_coefficients:
                        if coefficient not in [0.0, 1.0]:
                            print("*** WARNING ***\n")
                            print(
                                "Program work with uncontractred bases set only."
                            )
                            stop_verification_contraction = True

                line_gto += primitive_amount
                total_primitives += (
                    primitive_amount * l[spatial_primitive][primitive_type]
                )

                # To pass the next element
                if datafile[number_line + line_gto] == "\n":
                    line_gto += 1
                    stop_read_exponents = True

                if (
                    primitive_amount > 1
                    and datafile[number_line + line_gto].split()[0]
                    == primitive_type
                ):
                    # Skip following lines until new primitive type
                    # primitive_amount - 1: Number of blocks to skip
                    # primitive_amount + 1: Number of lines by block
                    line_gto += (primitive_amount - 1) * (primitive_amount + 1)
                    # To pass the next element, if there isn't another primitive type
                    if datafile[number_line + line_gto] == "\n":
                        line_gto += 1
                        stop_read_exponents = True
                elif (
                    primitive_amount == 1
                    and datafile[number_line + line_gto].split()[0]
                    == primitive_type
                ):
                    # TODO when the exponents for the same primitive type
                    # TODO are separated by a one line with same type and amount primitive
                    line_gto += 1
                elif (
                    primitive_amount == 1
                    and datafile[number_line + line_gto].split()[0]
                    != primitive_type
                ):
                    pass

                # primitive.append(
                #     [
                #         primitive_type,
                #         primitive_amount,
                #         exponents,
                #     ]
                # )
                primitive[primitive_type] = exponents

            t_a_exp.append(primitive)
            # Verication: there is another element
            if len(datafile[number_line + line_gto].split()) != 2:
                stop_gto = True
            else:
                line_gto += 1

        while "[MO]" in element and stop_mo == False:
            if not stop_gto:
                raise ValueError("[GTO] section must be before section [MO]")
            # MO information
            #
            # [MO]
            # Sym = A
            # Ene = ...
            # Spin = Alpha
            # Occup = ...
            # MO Coefficients
            #
            sym = datafile[number_line + count_mo].split()[1]
            count_mo += 1
            e = convert_to_float("energy",
                                datafile[number_line + count_mo].split()[1])
            count_mo += 1
            spin = datafile[number_line + count_mo].split()[1]
            count_mo += 1
            occupation = convert_to_float("occupation",
                                    datafile[number_line + count_mo].split()[1])
            count_mo += 1

            mo_coefficients = [0.0] * total_primitives
            count_different_0 = 0
            for i in range(total_primitives):
                if (
                    datafile[number_line + count_mo + i].split()[0][0]
                    not in end_mo
                ):
                    mo_coefficients[
                        int(datafile[number_line + count_mo + i].split()[0])
                        - 1
                    ] = convert_to_float("coefficient",
                                        datafile[number_line + count_mo + i].split()[1])
                    count_different_0 += 1
                else:
                    break

            count_mo += total_primitives - (
                total_primitives - count_different_0
            )

            mo_information = {
                "energy": e,
                "spin": spin,
                "occupation": occupation,
                "coefficients": mo_coefficients,
            }
            mo.append(mo_information)

            num_mo += 1
            if num_mo > total_primitives:
                stop_mo = True
            elif "Sym=" != datafile[number_line + count_mo].split()[0]:
                raise ValueError("MO Coefficients must be a square matrix")

    content.close()

    # NOTE: Re--organize of coefficient when there are angular moments higher than p
    #       because the coefficients aren't ordered by ml, when is used the DALTON
    angular_moments =[il for dicti in t_a_exp for il, values in dicti.items() for i in values for j in range(l[spatial_primitive][il])]
    if "d" in angular_moments and spatial_primitive == "spherical":
        drv_scratch.write_output("*** There are l higher p, then is neccesary re-organize mo coefficient ***")
        drv_scratch.write_output("*** when is used DALTON, which will done ***")
        new_mo = []
        for imo in mo:
            coefficients = [0.0 for i in range(total_primitives)]
            count = 0
            ml = 0
            for iterator, coeff in enumerate(imo["coefficients"]):
                nml = l["spherical"][angular_moments[iterator]]
                if nml > 3:
                    coefficients[count + int(nml/2) + int((ml+1)/2)*np.power(-1,ml+1)] = coeff
                    ml += 1
                else:
                    coefficients[count] = coeff
                    count += 1
                if  ml == nml:
                    ml = 0
                    count += nml
            new_mo.append({"energy": imo["energy"], "spin": imo["spin"], "occupation": imo["occupation"], "coefficients": coefficients})
        mo = new_mo
    elif "d" in angular_moments and spatial_primitive != "spherical":
        drv_scratch.write_output("*** There are l higher p, then is neccesary re-organize mo coefficient ***")
        drv_scratch.write_output("when is used DALTON. Also, it is used cartessian primitive and in this case\
                the primitives isn't reorganized")

    # Encapsule t_a_exp like l_i_q_xyz
    t_a_exp_temp = t_a_exp
    t_a_exp = []
    count = 0
    for mol in l_i_q_xyz:
        atomic_basis = []
        for atom in mol:
            atomic_basis.append(t_a_exp_temp[count])
            count += 1
        t_a_exp.append(atomic_basis)

    # NOTE: It's neccesary to verified if primitives are cartessian or
    #       spherical because this effect split of the atomic orbitals.
    #       In the calculation integrals always calculate in cartessian
    #       form, then integrals is convert to spherical.
    # WARNING: When the basis set used is spherical then is neccesary
    #           recalculated the molecular orbital energies, because
    #           if not there will be errors in the calculations where
    #           is neccesary

    if spatial_primitive == "spherical":
        return l_i_q_xyz, t_a_exp, mo, False
    else:
        return l_i_q_xyz, t_a_exp, mo, True

if __name__ == "__main__":
    """
    Read .molden file
    """

    l_i_q_xyz, t_a_exp, mo, type_primitives = read_molden("LiH.molden")

    # print("\n primitive information \n",t_a_exp)
    # print("\n mo coefficientes \n",mo[0]['coefficients'])

