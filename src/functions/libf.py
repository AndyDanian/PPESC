import os
import sys
from pathlib import Path
from time import time
from typing import List, Set, Dict, Tuple, Optional

# Addres when execute from fock
FUNCTIONS_PATH = Path.cwd()

PARENT_PATH = FUNCTIONS_PATH.parent

sys.path.append(os.fspath(PARENT_PATH))
sys.path.append(os.fspath(PARENT_PATH / ("include")))  # This is neccesary by e_integral
sys.path.append(os.fspath(PARENT_PATH / ("io")))  # This is neccesary by wave_function

import numpy as np

# project PATH
from wave_function import *

# include
from response_parameters import *
from integrals_parameters import *
