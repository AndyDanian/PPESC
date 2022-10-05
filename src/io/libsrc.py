import os
import sys
from pathlib import Path
from time import time

import numpy as np

# Add folders to system path
PROJECT_DIR = Path.cwd().parent
sys.path.append(
    os.fspath(PROJECT_DIR / "include")
    )
sys.path.append(
    os.fspath(PROJECT_DIR / "functions")
    )
# sys.path.append(
#     os.fspath(PROJECT_DIR / "fock")
#     )


from atom import *
from scratch import *
from print_matrix import *

#include
from atomic_symbol import *
from integrals_parameters import *
from quantum_numbers import *

# io
from molden import *

# functions
from convert_array import *
from drv_time import *
