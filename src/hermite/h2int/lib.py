import os
import sys
from pathlib import Path
from time import time

import numpy as np

H1I_PATH = Path.cwd()

sys.path.append(
    os.fspath(H1I_PATH.parent / 'normalization')
    )
sys.path.append(
    os.fspath(H1I_PATH.parent / 'eij')
    )
sys.path.append(
    os.fspath(H1I_PATH.parent / 'refg')
    )
sys.path.append(
    os.fspath(H1I_PATH.parent / 'gaussian_multiplication')
    )

# Normalization
from normalization import *
# Eij coefficients
from eij import *
# Refg
from refg import *
# Gaussian multiplication
from gaussian_multiplication import *
