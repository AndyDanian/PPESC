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
    os.fspath(PROJECT / ("functions"))
)

import numpy as np

# Hermite folder
from f90recursives import *

# functions
from string_informations import *
