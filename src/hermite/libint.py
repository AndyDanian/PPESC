import os
import sys
from pathlib import Path
from time import time

#Addres when execute from fock
HERMITE_PATH = Path.cwd()

PROJECT = HERMITE_PATH.parent

sys.path.append(
    os.fspath(PROJECT)
    )
sys.path.append(
    os.fspath(PROJECT / ("include")) #This is neccesary by e_integral
    )
sys.path.append(
    os.fspath(PROJECT / ("io")) #This is neccesary by wave_function
    )
sys.path.append(
    os.fspath(PROJECT / ("functions"))
    )
sys.path.append(
    os.fspath(PROJECT / ("include"))
    )
sys.path.append(
    os.fspath(HERMITE_PATH / ("h1int"))
    )
sys.path.append(
    os.fspath(HERMITE_PATH / ("h2int"))
    )

import numpy as np

from wave_function import *

#functions
from print_matrix import *

# functions
from integral_parameters import *
from convert_array import *
from print_matrix import *
from string_informations import *
from drv_time import *

#include
from constants_cto_gto import *
from integrals_parameters import *

#io
from molden import *

# h1int: One--Body hermite integrals
from overlap import *
from nucpot import *   # Nucleu potential
from kinetic import *  # Kinectic energy
from angmom import *   # Angular momentum
from sd import *       # Spin dipolar
from fc import *       # Fermi--contact
from darwin import *   # Darwin
from massvelo import * # Massvelo
from nelfld import *   # Nuclear electric field gradient
from diplen import *   # Dipole lenght
from dipvel import *   # Dipole velocity
from pso import *      # Paramagnetic spin-orbit
from nstcgo import *   # Diamagnetic nuclear shielding tensor
from dnske import *    # Kinetic-energy correction to the diamagnetic contribution to nuclear shielding
from psoke import *    # Kinetic-energy correction to the paramagnetic spin-orbit to nuclear shielding
from psooz import *    # Orbital-Zeeman correction to the paramagnetic spin-orbit to nuclear shielding
from ozke import *     # Calculates the kinetic energy correction to the orbital Zeeman operator
from laplacian import * # Calculate different double derivatives
from spinorbit import * # Calculate SpinOrbit integrals

# h2int: Two--Body hermite integrals
from e2pot import *