# from intherm1.overlap import *
import os
import sys
from pathlib import Path

PROJECT_DIR = Path.cwd()

sys.path.append(
    os.fspath(PROJECT_DIR / "hermite")
    )
sys.path.append(
    os.fspath(PROJECT_DIR / "hermite/h1int")
    )
sys.path.append(
    os.fspath(PROJECT_DIR / "hermite/h2int")
    )
sys.path.append(
    os.fspath(PROJECT_DIR / "include")
    )
sys.path.append(
    os.fspath(PROJECT_DIR / "io")
    )
sys.path.append(
    os.fspath(PROJECT_DIR / "functions")
    )

import numpy as np

#include
from dicts import *
from quantum_numbers import *

# io
from molden import *

# functions
from convert_array import *
from print_matrix import *

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

# h2int: Two--Body hermite integrals
from e2pot import *
