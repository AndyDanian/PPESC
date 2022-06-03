import os
import sys
from pathlib import Path
from time import time

#Addres when execute from fock
H1I_PATH = Path.cwd()

PARENT_PATH = H1I_PATH.parent

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

from e_integral import * 
from wave_function import * 

from print_matrix import *