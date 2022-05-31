
import numpy as np
from scipy.special import hyp1f1

from lib1h import *


def R(t, mu, nu, n, p, PKx, PKy, PKz, Rpc):
    T = p * Rpc * Rpc
    pot = 0.0
    if t == mu == nu == 0:
        pot += (
            np.power(-2.0 * p, n)
            * hyp1f1(n + 0.5, n + 1.5, -T)
            / (2.0 * n + 1.0)
        )
    elif t == mu == 0:
        if nu > 1:
            pot += (nu - 1) * R(t, mu, nu - 2, n + 1, p, PKx, PKy, PKz, Rpc)
        pot += PKz * R(t, mu, nu - 1, n + 1, p, PKx, PKy, PKz, Rpc)
    elif t == 0:
        if mu > 1:
            pot += (mu - 1) * R(t, mu - 2, nu, n + 1, p, PKx, PKy, PKz, Rpc)
        pot += PKy * R(t, mu - 1, nu, n + 1, p, PKx, PKy, PKz, Rpc)
    else:
        if t > 1:
            pot += (t - 1) * R(t - 2, mu, nu, n + 1, p, PKx, PKy, PKz, Rpc)
        pot += PKx * R(t - 1, mu, nu, n + 1, p, PKx, PKy, PKz, Rpc)
    return pot


def nuclear_attraction(
    i, k, m, j, l, n, e, f, g, alpha, beta, Ax, Ay, Az, Bx, By, Bz, Kx, Ky, Kz
):
    """
    Recurrence to calculate the integrate that include 1/r

    Equation 9.9.32 from Molecular Electronic-Structure Theory. T Helgaker, et al. 
    """
    p = alpha + beta

    Px = alpha * Ax + beta * Bx
    Px = Px / p
    Py = alpha * Ay + beta * By
    Py = Py / p
    Pz = alpha * Az + beta * Bz
    Pz = Pz / p

    Rpk = np.linalg.norm([Px - Kx, Py - Ky, Pz - Kz])

    suma = 0.0
    for t in range(i + j + 1):
        for mu in range(k + l + 1):
            for nu in range(m + n + 1):
                suma += (
                    E(i, j, t, Ax - Bx, alpha, beta)
                    * E(k, l, mu, Ay - By, alpha, beta)
                    * E(m, n, nu, Az - Bz, alpha, beta)
                    * R(
                        t + e,
                        mu + f,
                        nu + g,
                        0,
                        p,
                        Px - Kx,
                        Py - Ky,
                        Pz - Kz,
                        Rpk,
                    )
                )

    return suma * (-1) ** (e + f + g)
