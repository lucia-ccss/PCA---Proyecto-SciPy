# condiciones_contorno.py
"""Condiciones de contorno en los contactos Schottky (x=0) y óhmico (x=L)."""
import numpy as np

from equilibrio import calcular_Nc


def phi_contacto_ohmic(n0L, nisc, param):
    """Potencial en el contacto óhmico (referencia) [V]."""
    V_T = param.k_B * param.T
    return V_T * np.log(n0L / nisc)


def phi_contacto_schottky(nisc, V_app, param):
    """Potencial en el contacto Schottky [V]."""
    Nc  = calcular_Nc(param)
    V_T = param.k_B * param.T
    return V_T * np.log(Nc / nisc) - param.phi_B + V_app


def n_eq_schottky(param):
    """Densidad de equilibrio en el contacto Schottky [cm⁻³]."""
    Nc  = calcular_Nc(param)
    V_T = param.k_B * param.T
    return Nc * np.exp(-param.phi_B / V_T)


def bc_n_schottky(n1, phi0, phi1, n0eq, dx, param):
    """
    Densidad en j=0 después de actualizar los nodos interiores.
    Fórmula Robin con diferencias hacia adelante (notebook):

      n_0 = (v_e·n_e⁰(0) + D_e/dx·n_1) / (v_e + D_e/dx + μ_e/dx·(φ_1−φ_0))
    """
    num = param.v_e * n0eq + param.D_e / dx * n1
    den = param.v_e + param.D_e / dx + param.mu_e / dx * (phi1 - phi0)
    return num / den
