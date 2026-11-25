# main.py
"""
Simulación de la juntura GaAs–Au (Schottky–Óhmico) en 1D.
Resuelve las ecuaciones de continuidad y Poisson de forma acoplada
hasta alcanzar el estado estacionario para cada tensión aplicada.
Salida: curva I-V y perfiles de portadores / potencial.
"""
from pathlib import Path

import numpy as np

from lectura_parametros import cargar_parametros
from configurar_parametros import configurar_parametros
from validacion_param import validar
from equilibrio import calcular_nisc, n_equilibrio
from condiciones_contorno import phi_contacto_ohmic, phi_contacto_schottky
from simulacion import simular_vapp, barrido_tension, perfiles_en_tensiones
from graficos import graficar_curva_IV, graficar_perfiles_y_IV, mostrar

_DIR = Path(__file__).resolve().parent


def main():
    # ── Parámetros ───────────────────────────────────────────────────
    param = cargar_parametros(_DIR / "constantes.txt", _DIR / "parametros.txt")
    configurar_parametros(param)      # permite cambiar T, L, N (CLI o interactivo)

    dx = param.L / param.N
    dt = 0.4 * dx**2 / param.D_e

    # ── Validación ───────────────────────────────────────────────────
    r = validar(param, dt)
    print(f"Parámetros válidos.  D·dt/dx² = {r:.4f}  (límite 0.5)\n")

    # ── Magnitudes de equilibrio ─────────────────────────────────────
    nisc  = calcular_nisc(param)
    n0_L  = n_equilibrio(param.C, nisc)
    phi_N = phi_contacto_ohmic(n0_L, nisc, param)

    print(f"n_isc  = {nisc:.3e} cm⁻³")
    print(f"n0_L   = {n0_L:.3e} cm⁻³")
    print(f"phi_N  = {phi_N:.4f} V")
    print(f"phi_0  = {phi_contacto_schottky(nisc, 0, param):.4f} V  (en equilibrio)\n")

    # ── Condición inicial: densidad de equilibrio uniforme ───────────
    n_init = np.full(param.N + 1, n0_L)

    # ── Equilibrio inicial (warm start base) ─────────────────────────
    print("Calculando equilibrio inicial...")
    n_eq, _, _, _ = simular_vapp(0.0, n_init, nisc, n0_L, phi_N, dx, dt, param)

    # ── Barrido de tensión (curva I-V) ───────────────────────────────
    V_vals = np.linspace(-1.0, 0.5, 31)
    J_vals, n_ultimo = barrido_tension(V_vals, n_eq, nisc, n0_L, phi_N, dx, dt, param)

    # ── Perfiles en algunas tensiones ────────────────────────────────
    V_perfiles = [0.0, 0.3, -0.5]
    perfiles = perfiles_en_tensiones(V_perfiles, n_ultimo, nisc, n0_L, phi_N,
                                     dx, dt, param)

    # ── Gráficos ─────────────────────────────────────────────────────
    x_um = np.linspace(0, param.L * 1e4, param.N + 1)
    graficar_curva_IV(V_vals, J_vals, archivo=_DIR / "curva_IV.png")
    graficar_perfiles_y_IV(x_um, V_perfiles, perfiles, V_vals, J_vals,
                           archivo=_DIR / "perfiles_y_IV.png")
    mostrar()


if __name__ == "__main__":
    main()
