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

import numpy as np

# Normalization
from normalization import *
# Eij coefficients
from eij import *
# R1efg
from nuclear_attraction import *
# R2efg
from electron_repulsion import *
#functions
from string_informations import *
