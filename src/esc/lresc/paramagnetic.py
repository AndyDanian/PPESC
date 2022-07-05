from libl import *


def run_shielding_paramagnetic(wf: wave_function = None, paramagnetic: list = None,
                                atom: list = None,
                                principal_propagator_approximation: str = None,
                                driver_time: drv_time = None, verbose: int = 0):
    """
    Calculate all paramagnetic values
    Args:
        wf (wave_function): Wave function object
        paramagneitc (list): Paramagnetic amount to calculate
        atom (list): Atom index to calculate the shielding
        principal_propagator_approximation (str): Approximation
        driver_time (drv_time): Object to manage the time calculation
        verbose (int): Print level
    """
    start = time()

    for a in atom:
        gauge = wf.coordinates[a]
        fckin = ["fc " + str(1 + a) + ", kinetic"]

        sd1xddxx = ["sd " + str(1 + 3*a) + " 1, laplacian 1"] #In old LRESC used DPTOVL, which according of DALTON
        sd1xddxy = ["sd " + str(1 + 3*a) + " 1, laplacian 4"] #source is -1/4*ddi, !falta las cruzadas
        sd1xddxz = ["sd " + str(1 + 3*a) + " 1, laplacian 5"]
        sd1yddxx = ["sd " + str(1 + 3*a) + " 2, laplacian 1"]
        sd1yddxy = ["sd " + str(1 + 3*a) + " 2, laplacian 4"]
        sd1yddxz = ["sd " + str(1 + 3*a) + " 2, laplacian 5"]
        sd1zddxx = ["sd " + str(1 + 3*a) + " 3, laplacian 1"]
        sd1zddxy = ["sd " + str(1 + 3*a) + " 3, laplacian 4"]
        sd1zddxz = ["sd " + str(1 + 3*a) + " 3, laplacian 5"]

        sd2xddyx = ["sd " + str(2 + 3*a) + " 1, laplacian 4"] # dxy
        sd2xddyy = ["sd " + str(2 + 3*a) + " 1, laplacian 2"] # dyy
        sd2xddyz = ["sd " + str(2 + 3*a) + " 1, laplacian 6"] # dyz
        sd2yddyx = ["sd " + str(2 + 3*a) + " 2, laplacian 4"]
        sd2yddyy = ["sd " + str(2 + 3*a) + " 2, laplacian 2"]
        sd2yddyz = ["sd " + str(2 + 3*a) + " 2, laplacian 6"]
        sd2zddyx = ["sd " + str(2 + 3*a) + " 3, laplacian 4"]
        sd2zddyy = ["sd " + str(2 + 3*a) + " 3, laplacian 2"]
        sd2zddyz = ["sd " + str(2 + 3*a) + " 3, laplacian 6"]

        sd3xddzx = ["sd " + str(3 + 3*a) + " 1, laplacian 5"] # dxz
        sd3xddzy = ["sd " + str(3 + 3*a) + " 1, laplacian 6"] # dyz
        sd3xddzz = ["sd " + str(3 + 3*a) + " 1, laplacian 3"] # dzz
        sd3yddzx = ["sd " + str(3 + 3*a) + " 2, laplacian 5"]
        sd3yddzy = ["sd " + str(3 + 3*a) + " 2, laplacian 6"]
        sd3yddzz = ["sd " + str(3 + 3*a) + " 2, laplacian 3"]
        sd3zddzx = ["sd " + str(3 + 3*a) + " 3, laplacian 5"]
        sd3zddzy = ["sd " + str(3 + 3*a) + " 3, laplacian 6"]
        sd3zddzz = ["sd " + str(3 + 3*a) + " 3, laplacian 3"]

        angmomxpsoke1 = ["angmom 1, psoke " + str(1 + a*3)]
        angmomypsoke2 = ["angmom 2, psoke " + str(2 + a*3)]
        angmomzpsoke3 = ["angmom 3, psoke " + str(3 + a*3)]
        psoxozke1 = ["pso " + str(1 + 3*a) + ", ozke 1"]
        psoyozke2 = ["pso " + str(2 + 3*a) + ", ozke 2"]
        psozozke3 = ["pso " + str(3 + 3*a) + ", ozke 3"]
        angmomxfcspinox = ["angmom 1, fc " + str(1 + a) + ", spinorbit 1"]
        angmomyfcspinoy = ["angmom 2, fc " + str(1 + a) + ", spinorbit 2"]
        angmomzfcspinoz = ["angmom 3, fc " + str(1 + a) + ", spinorbit 3"]

        angmomxsd1xspinox = ["angmom 1, sd " + str(1 + 3*a) + "1, spinorbit 1"]
        angmomysd1yspinoy = ["angmom 1, sd " + str(1 + 3*a) + "2, spinorbit 2"]
        angmomzsd1zspinoz = ["angmom 1, sd " + str(1 + 3*a) + "3, spinorbit 3"]
        angmomxsd2xspinox = ["angmom 2, sd " + str(2 + 3*a) + "1, spinorbit 1"]
        angmomysd2yspinoy = ["angmom 2, sd " + str(2 + 3*a) + "2, spinorbit 2"]
        angmomzsd2zspinoz = ["angmom 2, sd " + str(2 + 3*a) + "3, spinorbit 3"]
        angmomxsd3xspinox = ["angmom 3, sd " + str(3 + 3*a) + "1, spinorbit 1"]
        angmomysd3yspinoy = ["angmom 3, sd " + str(3 + 3*a) + "2, spinorbit 2"]
        angmomzsd3zspinoz = ["angmom 3, sd " + str(3 + 3*a) + "3, spinorbit 3"]

        angmomxpso1massvelo = ["angmom 1, pso " + str(1 + 3*a) + ", massvelo"]
        angmomypso2massvelo = ["angmom 2, pso " + str(2 + 3*a) + ", massvelo"]
        angmomzpso3massvelo = ["angmom 3, pso " + str(3 + 3*a) + ", massvelo"]
        angmomxpso1darwin = ["angmom 1, pso " + str(1 + 3*a) + ", darwin"]
        angmomypso2darwin = ["angmom 2, pso " + str(2 + 3*a) + ", darwin"]
        angmomzpso3darwin = ["angmom 3, pso " + str(3 + 3*a) + ", darwin"]

    if verbose > 10:
        driver_time.add_name_delta_time(name = "Paramagnetic Calculations", delta_time = (time() - start))