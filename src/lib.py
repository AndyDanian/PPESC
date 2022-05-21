# from intherm1.overlap import *
import os
import sys
from pathlib import Path

PROJECT_DIR = Path.cwd()

print("PROJECT_DIR ",PROJECT_DIR)

sys.path.append(
    os.fspath(PROJECT_DIR / "h1int")
    )
sys.path.append(
    os.fspath(PROJECT_DIR / "io")
    )
sys.path.append(
    os.fspath(PROJECT_DIR / "functions")
    )

import numpy as np

# io
from molden import *
from print_matrix import *

# functions
from convert_array import *

# One--Body hermite integrals
from overlap import *
from pot import *
from kinetic import *
from angmom import *
