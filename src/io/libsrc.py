import os
import sys
from pathlib import Path
from time import time

PROJECT_DIR = Path.cwd().parent

# sys.path.append(
#     os.fspath(PROJECT_DIR / "hermite")
#     )
# sys.path.append(
#     os.fspath(PROJECT_DIR / "hermite/h1int")
#     )
# sys.path.append(
#     os.fspath(PROJECT_DIR / "hermite/h2int")
#     )
sys.path.append(
    os.fspath(PROJECT_DIR / "include")
    )
sys.path.append(
    os.fspath(PROJECT_DIR / "functions")
    )
sys.path.append(
    os.fspath(PROJECT_DIR / "fock")
    )

import numpy as np

from atom import *
from scratch import *
from print_matrix import *

#include
from integrals_parameters import *
from quantum_numbers import *
from atomic_symbol import *

# io
from molden import *

# functions
from convert_array import *
from drv_time import *

#fock
#from fock import *