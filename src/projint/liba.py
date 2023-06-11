import os
import sys
from pathlib import Path
from time import time
from typing import List, Set, Dict, Tuple, Optional

import numpy as np

# Addres when execute from fock
# This module used different integrals to build one, used:
#   <phi|A|Psi><Psi|B|phi>   phi: atomic orbital
#                            Psi: Molecular orbital
PROJINT_PATH = Path.cwd()

PROJECT = PROJINT_PATH.parent

sys.path.append(os.fspath(PROJECT))
# sys.path.append(os.fspath(PROJECT / ("include")))  # This is neccesary by e_integral
sys.path.append(os.fspath(PROJECT / ("io")))  # This is neccesary by wave_function
# sys.path.append(os.fspath(PROJECT / ("functions")))
sys.path.append(os.fspath(PROJECT / ("hermite")))
sys.path.append(os.fspath(PROJECT / ("hermite/h1int")))
sys.path.append(os.fspath(PROJECT / ("hermite/h2int")))

# project PATH
from wave_function import *

# hermite
from e_integral import *

# functions
from drv_time import *
