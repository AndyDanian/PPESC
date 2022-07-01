from libint import *


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
    verbose: int = 0,
    dalton_normalization: bool = None,
    driver_time: drv_time = None,
    atom: list = None,
    magnetic_xyz: int = None,
    spatial_sym: int = None,
    r_gauge: list = None,
    r_dipole: list = None,
):
    """
    Initialize the calculate of the integral choose

    Args:
        integrla (str): Integral name to calculate
        verbose (int): verbose level for integral calculation
        drive_time (drv_object): Object to manage the time

    """

    if name.lower() == "overlap":
        integral: list = overlap(coord, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time)

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
        integral: list = nucpot(charge, atom, coord, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time)
    elif name.lower() == "kinetic":
        integral: list = kinetic(coord, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time)
    elif name.lower() == "angmom":
        integral: list = angmom(coord, r_gauge, magnetic_xyz, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time)
    elif name.lower() == "sd":
        integral: list = sd(coord, magnetic_xyz, spatial_sym, atom, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time)
    elif name.lower() == "fc":
        integral: list = fc(coord, atom, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time)
    elif name.lower() == "darwin":
        integral: list = darwin(charge, coord, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time)
    elif name.lower() == "massvelo":
        integral: list = massvelo(coord, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time)
    elif name.lower() == "nelfld":
        integral: list = nelfld(coord, spatial_sym, atom, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time)
    elif name.lower() == "diplen":
        integral: list = diplen(coord, magnetic_xyz, r_dipole, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time)
    elif name.lower() == "dipvel":
        integral: list = dipvel(coord, magnetic_xyz, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time)
    elif name.lower() == "pso":
        integral: list = pso(coord, spatial_sym, atom, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time)
    elif name.lower() == "nstcgo":
        integral: list = nstcgo(coord, r_gauge, spatial_sym, magnetic_xyz, atom, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time)
    elif name.lower() == "dnske":
        integral: list = dnske(coord, r_gauge, spatial_sym, magnetic_xyz, atom, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time)
    elif name.lower() == "psoke":
        integral: list = psoke(coord, spatial_sym, atom, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time)
    elif name.lower() == "psooz":
        integral: list = psooz(coord, r_gauge, spatial_sym, magnetic_xyz, atom, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time)
    elif name.lower() == "ozke":
        integral: list = ozke(coord, r_gauge, magnetic_xyz, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time)
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
