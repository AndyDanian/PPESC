import os
import sys
from pathlib import Path
from time import time

#Addres when execute from hxint
H1I_PATH = Path.cwd()

HERMITE_PATH = H1I_PATH.parent
PROJECT = HERMITE_PATH.parent

sys.path.append(
    os.fspath(HERMITE_PATH)
    )
#
sys.path.append(
    os.fspath(HERMITE_PATH / ("h1int"))
    )
sys.path.append(
    os.fspath(HERMITE_PATH / ("h2int"))
    )

import numpy as np

# Normalization
from normalization import *
# Eij coefficients
from eij import *
# Refg
from nuclear_attraction import *
# Gaussian multiplication
from gaussian_multiplication import *
# functions
from string_informations import *
