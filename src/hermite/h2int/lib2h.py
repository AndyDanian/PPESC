import os
import sys
from pathlib import Path
from time import time

import numpy as np

# Addres when execute from hxint
H2I_PATH = Path.cwd()

HERMITE_PATH = H2I_PATH.parent
PROJECT = HERMITE_PATH.parent

sys.path.append(os.fspath(HERMITE_PATH))
#
sys.path.append(os.fspath(PROJECT / ("io")))
# These modules are neccesaries to run e2pot into h2int
sys.path.append(os.fspath(PROJECT / ("functions")))
sys.path.append(os.fspath(PROJECT / ("include")))
# Hermite folder
from f90recursives import *

# io
from wave_function import *
