#codigo que revisa que los parametros del material tengan sentido físico

# validacion_param.py
"""Chequeos de consistencia de parámetros antes de correr la simulación."""
import warnings
import numpy as np

K_B = 8.617e-5  # constante de Boltzmann [eV/K]


def _positivo(nombre, valor):
    if not np.isfinite(valor) or valor <= 0:
        raise ValueError(f"{nombre} debe ser finito y > 0 (recibido {valor!r})")


def validar_fisicos(p):
    """p: módulo o objeto con los atributos de parametros.py."""
    for nombre in ("mu_e", "D_e", "tau_e", "epsilon_r", "E_g", "T", "L"):
        _positivo(nombre, getattr(p, nombre))
    if p.C < 0:
        raise ValueError(f"El dopaje C no puede ser negativo (recibido {p.C})")
    if not isinstance(p.N, (int, np.integer)) or p.N < 3:
        raise ValueError(f"N debe ser un entero >= 3 (recibido {p.N!r})")


def validar_einstein(mu, D, T, rtol=0.1):
    """Relación de Einstein: D = mu * kT/q. Avisa si D y mu son inconsistentes."""
    D_esperado = mu * K_B * T
    err = abs(D - D_esperado) / D_esperado
    if err > rtol:
        warnings.warn(f"D_e={D} difiere {err:.1%} de mu*kT/q={D_esperado:.1f} cm²/s")
    return err

# revisar esta función
def validar_estabilidad(D, dt, dx, mu=None, E_max=None):
    """Estabilidad del esquema FTCS de continuidad.py."""
    r = D * dt / dx**2
    if r > 0.5:
        raise ValueError(
            f"Inestable: D*dt/dx² = {r:.3g} > 1/2. "
            f"Usá dt <= {0.5 * dx**2 / D:.3e} s"
        )
    # Péclet de malla: con diferencias centradas, si es > 1 aparecen oscilaciones
    if mu is not None and E_max is not None:
        Pe = mu * abs(E_max) * dx / (2 * D)
        if Pe > 1:
            warnings.warn(f"Péclet de malla = {Pe:.3g} > 1: posibles oscilaciones; refiná dx")
    return r


def validar_arrays(N, **arrays):
    """Uso: validar_arrays(N, n=n, phi=phi, G=G)"""
    for nombre, a in arrays.items():
        a = np.asarray(a)
        if a.shape != (N + 1,):
            raise ValueError(f"{nombre} tiene forma {a.shape}, se esperaba ({N + 1},)")
        if not np.all(np.isfinite(a)):
            raise ValueError(f"{nombre} contiene NaN o inf")


def validacion_de_parametros(p, dt):
    dx = p.L / p.N
    validar_fisicos(p)
    validar_einstein(p.mu_e, p.D_e, p.T)
    _positivo("dt", dt)
    return validar_estabilidad(p.D_e, dt, dx)
