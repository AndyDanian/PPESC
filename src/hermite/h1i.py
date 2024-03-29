from libint import *


def h1i(
    name: str,
    driver_time: drv_time,
    charge: list[float] = [],
    coord: list[list[float]] = [],
    exp: list[float] = [],
    center: list[int] = [],
    lx: list[int] = [],
    ly: list[int] = [],
    lz: list[int] = [],
    atom: list[int] = [],
    magnetic_xyz: int = -1,
    spatial_sym: int = -1,
    r_gauge: list[float] = [],
    r_dipole: list[float] = [],
    dalton_normalization: bool = False,
    verbose: int = 0,
):
    """
    Initialize the calculate of the integral choose

    Args:
        name (str): Integral name to calculate
        drive_time (drv_object): Object to manage the time
        charge (list): Atomic charge
        coord (list): Atomic Coordinate
        exp (list): Exponential of the primitives
        center (list): Where is centering of primitives
        lx, ly, lz (list): Component spatial of angular momentum component
        magnetic_xyz (int): Component spatial of the magnetic field
        spatial_sym (int): Spatial symmetry
        r_gauge (list): Gauge origen
        r_dipole (list): Dipole origen
        verbose (int): verbose level for integral calculation
    """

    if name.lower() == "overlap":
        integral: list = overlap(
            coord, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time
        )

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
        integral = nucpot(
            charge,
            atom,
            coord,
            exp,
            center,
            lx,
            ly,
            lz,
            verbose,
            dalton_normalization,
            driver_time,
        )
    elif name.lower() == "kinetic":
        integral = kinetic(
            coord, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time
        )
    elif name.lower() == "angmom":
        integral = angmom(
            coord,
            r_gauge,
            magnetic_xyz,
            exp,
            center,
            lx,
            ly,
            lz,
            verbose,
            dalton_normalization,
            driver_time,
        )
    elif name.lower() == "sd":
        integral = sd(
            coord,
            magnetic_xyz,
            spatial_sym,
            atom,
            exp,
            center,
            lx,
            ly,
            lz,
            verbose,
            dalton_normalization,
            driver_time,
        )
    elif name.lower() == "fc":
        integral = fc(
            coord,
            atom,
            exp,
            center,
            lx,
            ly,
            lz,
            verbose,
            dalton_normalization,
            driver_time,
        )
    elif name.lower() == "darwin":
        integral = darwin(
            charge,
            coord,
            exp,
            center,
            lx,
            ly,
            lz,
            verbose,
            dalton_normalization,
            driver_time,
        )
    elif name.lower() == "massvelo":
        integral = massvelo(
            coord, exp, center, lx, ly, lz, verbose, dalton_normalization, driver_time
        )
    elif name.lower() == "nelfld":
        integral = nelfld(
            coord,
            spatial_sym,
            atom,
            exp,
            center,
            lx,
            ly,
            lz,
            verbose,
            dalton_normalization,
            driver_time,
        )
    elif name.lower() == "diplen":
        integral = diplen(
            coord,
            magnetic_xyz,
            r_dipole,
            exp,
            center,
            lx,
            ly,
            lz,
            verbose,
            dalton_normalization,
            driver_time,
        )
    elif name.lower() == "dipvel":
        integral = dipvel(
            coord,
            magnetic_xyz,
            exp,
            center,
            lx,
            ly,
            lz,
            verbose,
            dalton_normalization,
            driver_time,
        )
    elif name.lower() == "pso":
        integral = pso(
            coord,
            spatial_sym,
            atom,
            exp,
            center,
            lx,
            ly,
            lz,
            verbose,
            dalton_normalization,
            driver_time,
        )
    elif name.lower() == "nstcgo":
        integral = nstcgo(
            coord,
            r_gauge,
            spatial_sym,
            magnetic_xyz,
            atom,
            exp,
            center,
            lx,
            ly,
            lz,
            verbose,
            dalton_normalization,
            driver_time,
        )
    elif name.lower() == "dnske":
        integral = dnske(
            coord,
            r_gauge,
            spatial_sym,
            magnetic_xyz,
            atom,
            exp,
            center,
            lx,
            ly,
            lz,
            verbose,
            dalton_normalization,
            driver_time,
        )
    elif name.lower() == "psoke":
        integral = psoke(
            coord,
            spatial_sym,
            atom,
            exp,
            center,
            lx,
            ly,
            lz,
            verbose,
            dalton_normalization,
            driver_time,
        )
    elif name.lower() == "psooz":
        integral = psooz(
            coord,
            r_gauge,
            spatial_sym,
            magnetic_xyz,
            atom,
            exp,
            center,
            lx,
            ly,
            lz,
            verbose,
            dalton_normalization,
            driver_time,
        )
    elif name.lower() == "ozke":
        integral = ozke(
            coord,
            r_gauge,
            magnetic_xyz,
            exp,
            center,
            lx,
            ly,
            lz,
            verbose,
            dalton_normalization,
            driver_time,
        )
    elif name.lower() == "laplacian":
        integral = laplacian(
            coord,
            magnetic_xyz,
            exp,
            center,
            lx,
            ly,
            lz,
            verbose,
            dalton_normalization,
            driver_time,
        )
    elif name.lower() == "pnstcgop":
        integral = pnstcgop(
            coord,
            r_gauge,
            spatial_sym,
            magnetic_xyz,
            atom,
            exp,
            center,
            lx,
            ly,
            lz,
            verbose,
            dalton_normalization,
            driver_time,
        )
    return integral
