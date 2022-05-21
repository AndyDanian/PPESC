import linecache
import numpy as np


def validate_file(file_molden):
    """
    Validate the file.
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
                print("Molden's file is valid.")
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


def read_molden(file_molden, verbose=21):
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
    # Verification of the file
    validate_file(file_molden)

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
        "cartessian": {"s": 1, "p": 3, "d": 6, "f": 10, "g": 15},
        "spherical": {"s": 1, "p": 3, "d": 5, "f": 7, "g": 9},
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
        while "[Atoms]" in element and stop_atoms == False:
            # Read atomic type and coordinates
            #
            # [Atoms] Au/Angstrom
            # Atom Symbol  Atomic Number  Charge  x  y  z
            # Atom Symbol  Atomic Number  Charge  x  y  z
            #
            if datafile[number_line].split()[1] != "AU":
                A_Bohr = 1.0 / 0.52917721
            if (
                amount_strings
                != len(datafile[number_line + 1 + numbers_atoms].split())
                and numbers_atoms > 0
            ):
                stop_atoms = True
                if len(l_i_q_xyz) == numbers_atoms:
                    print(
                        "*** Warning\n\
                    Each atom is a molecule. You can group atoms put on the same index\n\
                    into molden file in section [Atoms], second column. Example: \n\n\
                    O 1 8 0.0 0.0 0.0\n\
                    H 1 1 0.0 0.0 0.0\n\
                    H 1 1 0.0 0.0 0.0\n"
                    )
            else:
                molecule_information = [
                    datafile[number_line + 1 + numbers_atoms].split()[0]
                    + "  "
                    + datafile[number_line + 1 + numbers_atoms].split()[2]
                    + "  "
                    + str(
                        float(
                            datafile[number_line + 1 + numbers_atoms].split()[
                                3
                            ]
                        )
                        * A_Bohr
                    )
                    + "  "
                    + str(
                        float(
                            datafile[number_line + 1 + numbers_atoms].split()[
                                4
                            ]
                        )
                        * A_Bohr
                    )
                    + "  "
                    + str(
                        float(
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
                    line_gto += (primitive_amount - 1) * 7
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
            e = float(datafile[number_line + count_mo].split()[1])
            count_mo += 1
            spin = datafile[number_line + count_mo].split()[1]
            count_mo += 1
            occupation = float(datafile[number_line + count_mo].split()[1])
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
                    ] = float(datafile[number_line + count_mo + i].split()[1])
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

    if spatial_primitive == "spherical":
        return l_i_q_xyz, t_a_exp, mo, False
    else:
        return l_i_q_xyz, t_a_exp, mo, True