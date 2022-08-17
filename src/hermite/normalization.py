from numpy import sqrt, pi
import numpy as np

############################## Normalization factor ##############################
# N = (2*alpha/pi)^3/4*sqrt({[8*alpha]^(i+j+k) i!j!k!}/{[2i]![2j]![2k]!})

# https://numba.pydata.org/numba-doc/latest/developer/inlining.html?highlight=factorial
# Conflict with numba, produce that woulb be implemeted factorial function
def factorial(n):
    if n <= 0:
        return 1
    return n * factorial(n - 1)

def normalization(i: int = None, j: int = None, k: int = None, alpha: float = None, dalton_normalization: bool = False):
    if dalton_normalization:
        # Only with cto primitives, however with gto primitives use correct normalization
        return ((2.0 * alpha / pi) ** (3.0 / 4.0) * sqrt(( 4.0 * alpha ) ** (i + j + k)))
    else:
        return ((2.0 * alpha / pi) ** (3.0 / 4.0) *
        sqrt((( 8.0 * alpha ) ** (i + j + k) * factorial(i) * factorial(j) * factorial(k))
        / (factorial(2*i) * factorial(2*j) * factorial(2*k))))
