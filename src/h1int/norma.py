from math import factorial
from numpy import sqrt, exp, pi
from numpy import math


############################## Normalization factor ##############################
# N = (2*alpha/pi)^3/4*sqrt({[8*alpha]^(i+j+k) i!j!k!}/{[2i]![2j]![2k]!})

def normalization(i: int = None, j: int = None, k: int = None, alpha: float = None, dalton_normalization: bool = False): 
    if dalton_normalization:
        return ((2.0 * alpha / pi) ** (3.0 / 4.0) * sqrt(( 4.0 * alpha ) ** (i + j + k)))
    else:
        return ((2.0 * alpha / pi) ** (3.0 / 4.0) *
        sqrt((( 8.0 * alpha ) ** (i + j + k) * math.factorial(i) * math.factorial(j) * math.factorial(k))
        / (math.factorial(2*i) * math.factorial(2*j) * math.factorial(2*k))))



############################## Primitive Functions ##############################
# Phi -> S = SxSySy
Sx = lambda x, alpha, Rx: exp(-alpha * (x - Rx) ** 2)
Sy = lambda y, alpha, Ry: exp(-alpha * (y - Ry) ** 2)
Sz = lambda z, alpha, Rz: exp(-alpha * (z - Rz) ** 2)

# Phi -> P = PxPyPz
Px = lambda x, alpha, Rx: ((x - Rx) * Sx(x, alpha, Rx))
Py = lambda y, alpha, Ry: ((y - Ry) * Sy(y, alpha, Ry))
Pz = lambda z, alpha, Rz: ((z - Rz) * Sz(z, alpha, Rz))

# Phi -> D = DxDyDz
Dxx = lambda x, alpha, Rx: (Sx(x, alpha, Rx) * (x - Rx) ** 2)
Dyy = lambda y, alpha, Ry: (Sy(y, alpha, Ry) * (y - Ry) ** 2)
Dzz = lambda z, alpha, Rz: (Sz(z, alpha, Rz) * (z - Rz) ** 2)
