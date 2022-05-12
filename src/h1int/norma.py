from numpy import sqrt, exp, pi


############################## Normalization factor ##############################
# N = (2*alpha/pi)^3/4*sqrt({[8*alpha]^(i+j+k) i!j!k!}/{[2i]![2j]![2k]!})
NS = lambda alpha: (2.0 * alpha / pi) ** (3.0 / 4.0)
NP = lambda alpha: 2.0 * NS(alpha) * sqrt(alpha)
Norm = {0: NS, 1: NP}


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
