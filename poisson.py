#poisson
import numpy as np
from scipy.linalg import solve_banded


def resolver_poisson(n, p, phi0, phiN, dx, param):
    """
    Resuelve Poisson en los nodos interiores.

    Recibe:
        n, p   : arrays de largo N+1, densidades de electrones y huecos [cm⁻³]
        phi0   : potencial en x=0 [V]
        phiN   : potencial en x=L [V]
        dx     : paso espacial [cm]
        param  : objeto con los parámetros leídos en main.py
                 (usa param.C, param.e, param.epsilon_0, param.epsilon_r)

    Devuelve:
        phi    : array de largo N+1 con el potencial en toda la grilla [V]
    """
    N = len(n) - 1
    C = param.C
    if np.ndim(C) == 0:
        C = np.full(n.shape, C)
    b = dx**2 * param.e / (param.epsilon_0 * param.epsilon_r) * (p[1:-1] - n[1:-1] + C[1:-1])
    b[0]  += phi0          # contorno Schottky
    b[-1] += phiN          # contorno óhmico

    ab = np.zeros((3, N - 1))   # formato de bandas de solve_banded
    ab[0, 1:]  = -1             # superdiagonal
    ab[1, :]   =  2             # diagonal
    ab[2, :-1] = -1             # subdiagonal

    phi = np.empty(N + 1)
    phi[0], phi[-1] = phi0, phiN
    phi[1:-1] = solve_banded((1, 1), ab, b)
    return phi