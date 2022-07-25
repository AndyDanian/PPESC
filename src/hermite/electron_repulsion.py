import numpy as  np
from numba import njit

from lib2h import *

@njit
def electron_repulsion(
    # Alpha and Beta centers
    i, k, m, j, l, n, alpha, beta, Ax, Ay, Az, Bx, By, Bz,
    # Gamma and Delta centers
    u, v, w, x, y, z, gamma, delta, Cx, Cy, Cz, Dx, Dy, Dz
):
    """
    Recurrence to calculate the integrates the electron repulsion
                int phi_i phi_j 1/r phi_k phi_l dt

    Equation 9.9.33 from Molecular Electronic-Structure Theory. T Helgaker, et al.
    """
    p = alpha + beta  # center one gaussian, composite p's exponents (alpha, beta)
    q = gamma + delta # center one gaussian, composite q's exponents (gamma, delta)

    Px = alpha * Ax + beta * Bx
    Px = Px / p
    Py = alpha * Ay + beta * By
    Py = Py / p
    Pz = alpha * Az + beta * Bz
    Pz = Pz / p

    Qx = gamma * Cx + delta * Dx
    Qx = Qx / q
    Qy = gamma * Cy + delta * Dy
    Qy = Qy / q
    Qz = gamma * Cz + delta * Dz
    Qz = Qz / q

    pq = p*q/(p+q)

    RPQ = np.linalg.norm(np.array([Px - Qx, Py - Qy, Pz - Qz]))

    suma = 0.0
    for t in range(i + j + 1):
        for mu in range(k + l + 1):
            for nu in range(m + n + 1):
                for phi in range(u + x + 1):
                    for tau in range(v + y + 1):
                        for theta in range(w + z + 1):
                            suma += (
                                E(i, j, t, Ax - Bx, alpha, beta)
                                * E(k, l, mu, Ay - By, alpha, beta)
                                * E(m, n, nu, Az - Bz, alpha, beta)
                                * E(u, x, phi, Cx - Dx, gamma, delta)
                                * E(v, y, tau, Cy - Dy, gamma, delta)
                                * E(w, z, theta, Cz - Dz, gamma, delta)
                                * np.power(-1,phi+tau+theta)
                                * R(
                                    t + phi,
                                    mu + tau,
                                    nu + theta,
                                    0,
                                    pq,
                                    Px - Qx,
                                    Py - Qy,
                                    Pz - Qz,
                                    RPQ,
                                )
                            )
    suma *= 2.0*np.power(np.pi,2.5)/(p*q*np.sqrt(p + q))
    return suma