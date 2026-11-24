# equilibrio.py
"""Magnitudes intrínsecas y de equilibrio del GaAs."""
import numpy as np


def calcular_nisc(param):
    """
    Densidad intrínseca de portadores del GaAs [cm⁻³].

    n_isc = 2·(2π k_B T/h²)^(3/2) · (m_eeff·m_heff)^(3/4) · exp(−E_g/2k_BT)
    """
    prefactor = 2 * (2 * np.pi * param.k_BJ * param.T / param.h**2)**1.5
    masa      = (param.m_eeff * param.m_heff)**0.75
    E_g_J     = param.E_g * param.e
    expon     = np.exp(-E_g_J / (2 * param.k_BJ * param.T))
    return prefactor * masa * expon * 1e-6             # m⁻³ → cm⁻³


def calcular_Nc(param):
    """Densidad efectiva de estados en la banda de conducción [cm⁻³]."""
    return 2 * (2 * np.pi * param.m_eeff * param.k_BJ * param.T / param.h**2)**1.5 * 1e-6


def n_equilibrio(C, nisc):
    """Densidad de equilibrio: C/2 + sqrt((C/2)² + n_isc²) [cm⁻³]."""
    return C / 2 + np.sqrt((C / 2)**2 + nisc**2)
