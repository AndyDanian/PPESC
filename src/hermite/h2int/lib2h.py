import os
import sys
from pathlib import Path
from time import time

#from typing import List, Set, Dict, Tuple, Optional

#Addres when execute from hxint
H2I_PATH = Path.cwd()

HERMITE_PATH = H2I_PATH.parent
PROJECT = HERMITE_PATH.parent

sys.path.append(
    os.fspath(HERMITE_PATH)
    )
#
sys.path.append(
    os.fspath(PROJECT /("functions"))
    )
sys.path.append(
    os.fspath(PROJECT /("io"))
    )
sys.path.append(
    os.fspath(PROJECT /("include"))
    )
sys.path.append(
    os.fspath(PROJECT /("fock"))
    )
sys.path.append(
    os.fspath(HERMITE_PATH /("h1int"))
    )

import numpy as np

# Hermite folder
from f90recursives import *

#functions
from string_informations import *
from drv_time import *
#functions
from wave_function import *
