from lib1h import *

def gaussian_mult(
    i, k, m, j, l, n, Ax, Ay, Az, Bx, By, Bz, alpha, beta, Kx, Ky, Kz
):
    # <phi|delta(rk)|phi> = <phi|delta(x-Kx,y-Ky,z-Kz)|phi> =
    # (xk-Ax)^i(yk-Ay)^k(zk-Az)^mexp(-alpha*(rk-A)^2) *
    # (xk-Bx)^j(yk-By)^l(zk-Bz)^nexp(-alpha*(rk-B)^2)

    gij = np.exp(-alpha * (Kx - Ax) ** 2) * np.exp(-beta * (Kx - Bx) ** 2)
    gij = np.power((Kx - Ax), i) * np.power((Kx - Bx), j) * gij

    gkl = np.exp(-alpha * (Ky - Ay) ** 2) * np.exp(-beta * (Ky - By) ** 2)
    gkl = np.power((Ky - Ay), k) * np.power((Ky - By), l) * gkl

    gmn = np.exp(-alpha * (Kz - Az) ** 2) * np.exp(-beta * (Kz - Bz) ** 2)
    gmn = np.power((Kz - Az), m) * np.power((Kz - Bz), n) * gmn

    return gij * gkl * gmn