import os
import sys
from pathlib import Path
from time import time
from typing import List, Set, Dict, Tuple, Optional

#Addres when execute from fock
FOCK_PATH = Path.cwd()

PARENT_PATH = FOCK_PATH.parent

sys.path.append(
    os.fspath(PARENT_PATH)
    )
sys.path.append(
    os.fspath(PARENT_PATH / ("include")) #This is neccesary by e_integral
    )
sys.path.append(
    os.fspath(PARENT_PATH / ("io")) #This is neccesary by wave_function
    )
sys.path.append(
    os.fspath(PARENT_PATH / ("functions"))
    )
sys.path.append(
    os.fspath(PARENT_PATH / ("hermite"))
    )
sys.path.append(
    os.fspath(PARENT_PATH / ("hermite/h1int"))
    )
sys.path.append(
    os.fspath(PARENT_PATH / ("hermite/h2int"))
    )

import numpy as np

# project PATH
from wave_function import *

# hermite
from e_integral import *

# functions
from print_matrix import *
from string_informations import *
from drv_time import *