import os
import sys
from pathlib import Path
from time import time

#Addres when execute from hxint
H1I_PATH = Path.cwd()

PARENT_PATH = H1I_PATH.parent

sys.path.append(
    os.fspath(PARENT_PATH)
    )
#

import numpy as np

# Normalization
from normalization import *
# Eij coefficients
from eij import *
# Refg
from refg import *
# Gaussian multiplication
from gaussian_multiplication import *
