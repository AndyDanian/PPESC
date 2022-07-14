import os
import sys
from pathlib import Path
from time import time
from typing import List, Set, Dict, Tuple, Optional

#Addres when execute from fock
LRESC_PATH = Path.cwd()

PROJECT_PATH = LRESC_PATH.parent.parent

sys.path.append(
    os.fspath(PROJECT_PATH)
    )
sys.path.append(
    os.fspath(PROJECT_PATH / ("include")) #This is neccesary by e_integral
    )
sys.path.append(
    os.fspath(PROJECT_PATH / ("io")) #This is neccesary by wave_function
    )
sys.path.append(
    os.fspath(PROJECT_PATH / ("functions"))
    )
sys.path.append(
    os.fspath(PROJECT_PATH / ("response"))
    )
sys.path.append(
    os.fspath(PROJECT_PATH / ("hermite"))
    )
sys.path.append(
    os.fspath(PROJECT_PATH / ("hermite/h1int"))
    )
sys.path.append(
    os.fspath(PROJECT_PATH / ("hermite/h2int"))
    )
sys.path.append(
    os.fspath(PROJECT_PATH / ("average"))
    )

import numpy as np

# project PATH
from wave_function import *

# include
from lresc_parameters import *
from response_parameters import *

# response
from response import *

# functions
from print_matrix import *
from integral_parameters import *
from string_informations import *
from drv_time import *

# average
from average import *