import os
import sys
from pathlib import Path
from time import time

from typing import List, Set, Dict, Tuple, Optional

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
# R1efg
from nuclear_attraction import *
# R2efg
from electron_repulsion import *
