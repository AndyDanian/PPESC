from numpy import sqrt, pi
from numpy import math


############################## Normalization factor ##############################
# N = (2*alpha/pi)^3/4*sqrt({[8*alpha]^(i+j+k) i!j!k!}/{[2i]![2j]![2k]!})

def normalization(i: int = None, j: int = None, k: int = None, alpha: float = None, dalton_normalization: bool = False):
    if dalton_normalization:
        # Only with cto primitives, however with gto primitives use correct normalization
        return ((2.0 * alpha / pi) ** (3.0 / 4.0) * sqrt(( 4.0 * alpha ) ** (i + j + k)))
    else:
        return ((2.0 * alpha / pi) ** (3.0 / 4.0) *
        sqrt((( 8.0 * alpha ) ** (i + j + k) * math.factorial(i) * math.factorial(j) * math.factorial(k))
        / (math.factorial(2*i) * math.factorial(2*j) * math.factorial(2*k))))
