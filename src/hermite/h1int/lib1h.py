import os
import sys
from pathlib import Path
from time import time

import numpy as np

# Addres when execute from hxint
H1I_PATH = Path.cwd()

HERMITE_PATH = H1I_PATH.parent
PROJECT = HERMITE_PATH.parent

sys.path.append(os.fspath(HERMITE_PATH))
sys.path.append(os.fspath(PROJECT / ("functions")))
sys.path.append(os.fspath(PROJECT / ("io")))
#
# Hermite folder
from f90recursives import *

# io
from scratch import scratch

# functions
from drv_time import drv_time
