from lib import *


# class h1i:
#     def __init__(self,charge,
#             coord,
#             exp,
#             center,
#             lx,
#             ly,
#             lz,):
#         """
#         Manages the hermite one--body electronic integrals calculations

#         Args:
#         """

#         if not wf:
#             raise ValueError(
#                 "*** Error \n\n There isn't information in the  wave function object."
#             )

#         self.wf = wf

#     ##################################################################
#     # METHODS
#     ##################################################################
#     def mo_coeff(self):
#         """ """
#         return [mo["coefficients"] for mo in self.wf["mos"]]

#     def mo_occ(self):
#         return [mo["occupation"] for mo in self.wf["mos"]]


def h1i(
    charge: list = None,
    coord: list = None,
    exp: list = None,
    center: list = None,
    lx: list = None,
    ly: list = None,
    lz: list = None,
    name: str = None,
    output: int = 0,
    atom: list = None,
    magnetic_xyz: int = None, 
    spatial_sym: int = None,
    gauge: list = None,
    rdipole: list = None,
):
    """
    Initialize the calculate of the integral choose

    Args:
        integrla (str): Integral name to calculate
        output (int): Output level for integral calculation
        0: Nothing
        1: Print calculation time
        11: Print the integrals
        atoms (int): atom number to Rk
    """

    if name.lower() == "overlap":
        integral: list = overlap(coord, exp, center, lx, ly, lz, output)

        # Comprobation of the basis set
        # mo_integral = self.molecular_matrix(integral, sym)
        # mo_occ = self.mo_occ()
        # ne = 0
        # for i, occ in enumerate(mo_occ):
        #     ne += occ * mo_integral[i][i]
        # if abs(ne - sum(mo_occ)) < 1e-5:
        #     print("The basis set hold the electron number ", ne)
        # else:
        #     print(
        #         "The basis set low precision in hold the electron number ",
        #         ne,
        #     )
    elif name.lower() == "nucpot":
        integral: list = nucpot(charge, atom, coord, exp, center, lx, ly, lz, output)
    elif name.lower() == "kinetic":
        integral: list = kinetic(coord, exp, center, lx, ly, lz, output)
    elif name.lower() == "angmom":
        integral: list = angmom(coord, gauge, magnetic_xyz, exp, center, lx, ly, lz, output)
    elif name.lower() == "sd":
        integral: list = sd(coord, magnetic_xyz, spatial_sym, atom, exp, center, lx, ly, lz, output)
    elif name.lower() == "fc":
        integral: list = fc(coord, atom, exp, center, lx, ly, lz, output)
    elif name.lower() == "darwin":
        integral: list = darwin(charge, coord, exp, center, lx, ly, lz, output)
    elif name.lower() == "massvelo":
        integral: list = massvelo(coord, exp, center, lx, ly, lz, output)
    elif name.lower() == "nelfld":
        integral: list = nelfld(coord, spatial_sym, atom, exp, center, lx, ly, lz, output)
    elif name.lower() == "diplen":
        integral: list = diplen(coord, magnetic_xyz, rdipole, exp, center, lx, ly, lz, output)
    elif name.lower() == "dipvel":
        integral: list = dipvel(coord, magnetic_xyz, exp, center, lx, ly, lz, output)
    elif name.lower() == "pso":
        integral: list = pso(coord, spatial_sym, atom, exp, center, lx, ly, lz, output)
    elif name.lower() == "nstcgo":
        integral: list = nstcgo(coord, gauge, spatial_sym, magnetic_xyz, atom, exp, center, lx, ly, lz, output)
    elif name.lower() == "dnske":
        integral: list = dnske(coord, gauge, spatial_sym, magnetic_xyz, atom, exp, center, lx, ly, lz, output)
    elif name.lower() == "psoke":
        integral: list = psoke(coord, spatial_sym, atom, exp, center, lx, ly, lz, output)
    elif name.lower() == "psooz":
        integral: list = psooz(coord, gauge, spatial_sym, magnetic_xyz, atom, exp, center, lx, ly, lz, output)
    return integral

    # def molecular_matrix(self, integral: list = None, sym: str = None):
    #     """
    #     Calculate the average value of the integral in the MO basis

    #     Args:
    #         integral (array): array 2d with atomic integrals

    #     Return:
    #         average (float): average value of the integral
    #     """

    #     mo_coefficients = self.mo_coeff()

    #     total_nprim = len(mo_coefficients[0][:])
    #     matriz_integral = vector_to_matrix(total_nprim, integral, sym)

    #     # Integral in molecular basis set
    #     mo_integral = list(
    #         np.matmul(
    #             np.array(mo_coefficients),
    #             np.matmul(
    #                 np.array(matriz_integral), np.array(mo_coefficients).T
    #             ),
    #         )
    #     )

    #     return mo_integral
