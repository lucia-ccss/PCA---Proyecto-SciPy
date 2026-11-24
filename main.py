# main.py
"""
Simulación de la juntura GaAs–Au (Schottky–Óhmico) en 1D.
Resuelve las ecuaciones de continuidad y Poisson de forma acoplada
hasta alcanzar el estado estacionario para cada tensión aplicada.
Salida: curva I-V y perfiles de portadores / potencial.
"""

import ast
import operator
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import matplotlib.pyplot as plt

from validacion_param import validar, validar_arrays
from continuidad import paso_continuidad
from poisson import resolver_poisson


# ════════════════════════════════════════════════════════════════════
# Lectura de parámetros (parametros.txt) y constantes (constantes.txt)
# ════════════════════════════════════════════════════════════════════
# Formato de los .txt:  nombre = valor   # comentario
# El valor puede ser un número o una expresión aritmética simple
# (ej: 0.063 * 9.11e-31). Las líneas vacías y las que empiezan con
# '#' se ignoran.

_OPERADORES = {
    ast.Add: operator.add, ast.Sub: operator.sub,
    ast.Mult: operator.mul, ast.Div: operator.truediv,
    ast.Pow: operator.pow, ast.USub: operator.neg, ast.UAdd: operator.pos,
}


def _evaluar(nodo):
    """Evalúa de forma segura una expresión aritmética (sin usar eval)."""
    if isinstance(nodo, ast.Constant) and isinstance(nodo.value, (int, float)):
        return nodo.value
    if isinstance(nodo, ast.BinOp) and type(nodo.op) in _OPERADORES:
        return _OPERADORES[type(nodo.op)](_evaluar(nodo.left), _evaluar(nodo.right))
    if isinstance(nodo, ast.UnaryOp) and type(nodo.op) in _OPERADORES:
        return _OPERADORES[type(nodo.op)](_evaluar(nodo.operand))
    raise ValueError(f"Expresión no permitida: {ast.unparse(nodo)!r}")


def leer_txt(ruta):
    """Lee un archivo 'nombre = valor  # comentario' y devuelve un dict."""
    valores = {}
    with open(ruta, encoding="utf-8") as f:
        for num, linea in enumerate(f, start=1):
            linea = linea.split("#", 1)[0].strip()   # saco comentarios
            if not linea:
                continue
            if "=" not in linea:
                raise ValueError(f"{ruta}:{num}: falta '=' en {linea!r}")
            nombre, expr = (s.strip() for s in linea.split("=", 1))
            if not nombre.isidentifier():
                raise ValueError(f"{ruta}:{num}: nombre inválido {nombre!r}")
            valores[nombre] = _evaluar(ast.parse(expr, mode="eval").body)
    return valores


_DIR = Path(__file__).resolve().parent
_constantes = leer_txt(_DIR / "constantes.txt")
_material   = leer_txt(_DIR / "parametros.txt")

_repetidos = _constantes.keys() & _material.keys()
if _repetidos:
    raise ValueError(f"Parámetros definidos en ambos .txt: {sorted(_repetidos)}")

# Objeto con todos los parámetros como atributos (param.mu_e, param.T, ...).
# Se pasa como argumento a las funciones de continuidad.py, poisson.py y
# validacion_param.py.
param = SimpleNamespace(**_constantes, **_material)

# ════════════════════════════════════════════════════════════════════
# Magnitudes intrínsecas y de equilibrio
# ════════════════════════════════════════════════════════════════════

def calcular_nisc():
    """
    Densidad intrínseca de portadores del GaAs [cm⁻³].

    n_isc = 2·(2π k_B T/h²)^(3/2) · (m_eeff·m_heff)^(3/4) · exp(−E_g/2k_BT)
    """
    prefactor = 2 * (2 * np.pi * param.k_BJ * param.T / param.h**2)**1.5
    masa      = (param.m_eeff * param.m_heff)**0.75
    E_g_J     = param.E_g * param.e
    expon     = np.exp(-E_g_J / (2 * param.k_BJ * param.T))
    return prefactor * masa * expon * 1e-6             # m⁻³ → cm⁻³


def calcular_Nc():
    """Densidad efectiva de estados en la banda de conducción [cm⁻³]."""
    return 2 * (2 * np.pi * param.m_eeff * param.k_BJ * param.T / param.h**2)**1.5 * 1e-6


def n_equilibrio(C, nisc):
    """Densidad de equilibrio: C/2 + sqrt((C/2)² + n_isc²) [cm⁻³]."""
    return C / 2 + np.sqrt((C / 2)**2 + nisc**2)


# ════════════════════════════════════════════════════════════════════
# Condiciones de contorno
# ════════════════════════════════════════════════════════════════════

def phi_contacto_ohmic(n0L, nisc):
    """Potencial en el contacto óhmico (referencia) [V]."""
    V_T = param.k_B * param.T
    return V_T * np.log(n0L / nisc)


def phi_contacto_schottky(nisc, V_app):
    """Potencial en el contacto Schottky [V]."""
    Nc  = calcular_Nc()
    V_T = param.k_B * param.T
    return V_T * np.log(Nc / nisc) - param.phi_B + V_app


def n_eq_schottky(nisc):
    """Densidad de equilibrio en el contacto Schottky [cm⁻³]."""
    Nc  = calcular_Nc()
    V_T = param.k_B * param.T
    return Nc * np.exp(-phi_B / V_T)


def bc_n_schottky(n1, phi0, phi1, n0eq, dx):
    """
    Densidad en j=0 después de actualizar los nodos interiores.
    Fórmula Robin con diferencias hacia adelante (notebook):

      n_0 = (v_e·n_e⁰(0) + D_e/dx·n_1) / (v_e + D_e/dx + μ_e/dx·(φ_1−φ_0))
    """
    num = v_e * n0eq + param.D_e / dx * n1
    den = v_e + param.D_e / dx + param.mu_e / dx * (phi1 - phi0)
    return num / den


# ════════════════════════════════════════════════════════════════════
# Cálculo de corriente
# ════════════════════════════════════════════════════════════════════

def calcular_corriente(n, phi, dx):
    """Densidad de corriente promediada en x [A/cm²]. j_e = −e·J_e"""
    dphi  = (phi[1:] - phi[:-1]) / dx
    n_mid = 0.5 * (n[1:] + n[:-1])
    dn    = (n[1:] - n[:-1]) / dx
    J     = param.mu_e * dphi * n_mid - param.D_e * dn
    return -param.e * np.mean(J)


# ════════════════════════════════════════════════════════════════════
# Loop temporal: estado estacionario para un V_app dado
# ════════════════════════════════════════════════════════════════════

def simular_vapp(V_app, n_init, nisc, n0_L, phi_N, dx, dt,
                 G=0.0, tol=1e-6, max_iter=600_000):
    """
    Itera continuidad + Poisson hasta convergencia (estado estacionario).

    Retorna: n, phi, corriente [A/cm²], convergió (bool)
    """
    phi_0 = phi_contacto_schottky(nisc, V_app)
    n0_0  = n_eq_schottky(nisc)
    G_arr = np.broadcast_to(G, n_init.shape)

    n    = n_init.copy()
    phi  = resolver_poisson(n, np.zeros_like(n), phi_0, phi_N, dx, param)
    conv = False

    for _ in range(max_iter):
        n_prev = n.copy()

        # 1. avanzar nodos interiores
        n = paso_continuidad(n, phi, G_arr, dt, dx, v_e, n0_0, n0_L, param)

        # 2. reemplazar n[0] con la fórmula Robin del notebook
        n[0] = bc_n_schottky(n[1], phi_0, phi[1], n0_0, dx)

        # 3. resolver Poisson (p=0: minoritarios despreciables en GaAs n-type)
        phi = resolver_poisson(n, np.zeros_like(n), phi_0, phi_N, dx, param)

        # 4. criterio de convergencia
        delta = np.max(np.abs(n - n_prev)) / (np.max(n) + 1e-30)
        if delta < tol:
            conv = True
            break

    if not conv:
        print(f"  ⚠ No convergió en {max_iter} pasos (delta={delta:.2e})")

    return n, phi, calcular_corriente(n, phi, dx), conv


# ════════════════════════════════════════════════════════════════════
# Función principal
# ════════════════════════════════════════════════════════════════════

def main():
    dx = param.L / param.N
    dt = 0.4 * dx**2 / param.D_e

    # ── Validación ───────────────────────────────────────────────────
    r = validar(param, dt)
    print(f"Parámetros válidos.  D·dt/dx² = {r:.4f}  (límite 0.5)\n")

    # ── Magnitudes de equilibrio ─────────────────────────────────────
    nisc  = calcular_nisc()
    n0_L  = n_equilibrio(param.C, nisc)
    phi_N = phi_contacto_ohmic(n0_L, nisc)

    print(f"n_isc  = {nisc:.3e} cm⁻³")
    print(f"n0_L   = {n0_L:.3e} cm⁻³")
    print(f"phi_N  = {phi_N:.4f} V")
    print(f"phi_0  = {phi_contacto_schottky(nisc, 0):.4f} V  (en equilibrio)\n")

    # ── Condición inicial: densidad de equilibrio uniforme ───────────
    n_init = np.full(param.N + 1, n_equilibrio(param.C, nisc))

    # ── Equilibrio inicial (warm start base) ─────────────────────────
    print("Calculando equilibrio inicial...")
    n_actual, _, _, _ = simular_vapp(
        0.0, n_init, nisc, n0_L, phi_N, dx, dt
    )

    # ── Barrido de tensión ───────────────────────────────────────────
    V_vals = np.linspace(-1.0, 0.5, 31)
    J_vals = []                            # lista: se llena en el loop

    print(f"\n{'V_app':>8}  {'j [A/cm²]':>14}  {'conv':>5}")
    print("-" * 35)

    for V_app in V_vals:
        n_actual, phi_fin, j, conv = simular_vapp(
            V_app, n_actual, nisc, n0_L, phi_N, dx, dt   # warm start
        )
        J_vals.append(j)
        print(f"{V_app:+8.3f}  {j:+14.4e}  {'✓' if conv else 'X'}")

    J_vals = np.array(J_vals)             # convertir después del loop

    # ── Gráfica 1: curva I-V ─────────────────────────────────────────
    fig, ax = plt.subplots()
    ax.plot(V_vals, J_vals * 1e3, 'o-', lw=1.5, ms=4)
    ax.axhline(0, color='k', lw=0.5)
    ax.axvline(0, color='k', lw=0.5)
    ax.set_xlabel("$V_{App}$ [V]")
    ax.set_ylabel("$j_e$ [mA/cm²]")
    ax.set_title("Curva I–V: juntura GaAs–Au (Schottky–Óhmico)")
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    fig.savefig("curva_IV.png", dpi=150)
    print("\nGuardada: curva_IV.png")

    # ── Gráfica 2: perfiles en tres voltajes + I-V log ───────────────
    x_um = np.linspace(0, param.L * 1e4, param.N + 1)

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=False)
    etiquetas = [("$V_{App}=0$ V", 0.0),
                 ("$V_{App}=+0.3$ V", 0.3),
                 ("$V_{App}=-0.5$ V", -0.5)]

    n_base = n_actual.copy()
    for label, Vap in etiquetas:
        n_p, phi_p, _, _ = simular_vapp(
            Vap, n_base, nisc, n0_L, phi_N, dx, dt
        )
        axes[0, 0].semilogy(x_um, n_p, label=label)
        axes[1, 0].plot(x_um, phi_p, label=label)

    axes[0, 0].set_ylabel("$n_e$ [cm⁻³]")
    axes[0, 0].set_title("Densidad de electrones")
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].grid(True, alpha=0.4)
    axes[0, 0].set_xlabel("$x$ [µm]")

    axes[1, 0].set_ylabel("$\\phi$ [V]")
    axes[1, 0].set_title("Potencial electrostático")
    axes[1, 0].set_xlabel("$x$ [µm]")
    axes[1, 0].legend(fontsize=8)
    axes[1, 0].grid(True, alpha=0.4)

    axes[0, 1].plot(V_vals, J_vals * 1e3, 'o-', lw=1.5, ms=4, color='C1')
    axes[0, 1].axhline(0, color='k', lw=0.5)
    axes[0, 1].set_ylabel("$j_e$ [mA/cm²]")
    axes[0, 1].set_title("Curva I–V (lineal)")
    axes[0, 1].set_xlabel("$V_{App}$ [V]")
    axes[0, 1].grid(True, alpha=0.4)

    mask_pos = J_vals > 0
    mask_neg = J_vals < 0
    if mask_pos.any():
        axes[1, 1].semilogy(V_vals[mask_pos],  J_vals[mask_pos] * 1e3,
                            'o-', lw=1.5, ms=4, color='C1', label="$j > 0$")
    if mask_neg.any():
        axes[1, 1].semilogy(V_vals[mask_neg], -J_vals[mask_neg] * 1e3,
                            's--', lw=1.5, ms=4, color='C3', label="$|j|$ ($j<0$)")
    axes[1, 1].set_ylabel("|$j_e$| [mA/cm²]")
    axes[1, 1].set_title("Curva I–V (log)")
    axes[1, 1].set_xlabel("$V_{App}$ [V]")
    axes[1, 1].legend(fontsize=8)
    axes[1, 1].grid(True, alpha=0.4)

    fig.tight_layout()
    fig.savefig("perfiles_y_IV.png", dpi=150)
    print("Guardada: perfiles_y_IV.png")

    plt.show()


if __name__ == "__main__":
    main()