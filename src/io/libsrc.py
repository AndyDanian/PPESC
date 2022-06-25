# from intherm1.overlap import *
import os
import sys
from pathlib import Path

PROJECT_DIR = Path.cwd().parent

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
sys.path.append(
    os.fspath(PROJECT_DIR / "fock")
    )

import numpy as np

from atom import *

#include
from quantum_numbers import *
from atomic_symbol import *

# io
from molden import *

# functions
from convert_array import *
from print_matrix import *

#fock
from fock import *