# simulacion.py
"""
Loop temporal acoplado continuidad + Poisson hasta el estado estacionario,
y barrido en tensión aplicada para obtener la curva I-V.
"""
import numpy as np

from continuidad import paso_continuidad
from poisson import resolver_poisson
from condiciones_contorno import phi_contacto_schottky, n_eq_schottky, bc_n_schottky


def calcular_corriente(n, phi, dx, param):
    """Densidad de corriente promediada en x [A/cm²]. j_e = −e·J_e"""
    dphi  = (phi[1:] - phi[:-1]) / dx
    n_mid = 0.5 * (n[1:] + n[:-1])
    dn    = (n[1:] - n[:-1]) / dx
    J     = param.mu_e * dphi * n_mid - param.D_e * dn
    return -param.e * np.mean(J)


def simular_vapp(V_app, n_init, nisc, n0_L, phi_N, dx, dt, param,
                 G=0.0, tol=1e-6, max_iter=600_000):
    """
    Itera continuidad + Poisson hasta convergencia (estado estacionario).

    Retorna: n, phi, corriente [A/cm²], convergió (bool)
    """
    phi_0 = phi_contacto_schottky(nisc, V_app, param)
    n0_0  = n_eq_schottky(param)
    G_arr = np.broadcast_to(G, n_init.shape)

    n    = n_init.copy()
    phi  = resolver_poisson(n, np.zeros_like(n), phi_0, phi_N, dx, param)
    conv = False

    for _ in range(max_iter):
        n_prev = n.copy()

        # 1. avanzar nodos interiores
        n = paso_continuidad(n, phi, G_arr, dt, dx, param.v_e, n0_0, n0_L, param)

        # 2. reemplazar n[0] con la fórmula Robin del notebook
        n[0] = bc_n_schottky(n[1], phi_0, phi[1], n0_0, dx, param)

        # 3. resolver Poisson (p=0: minoritarios despreciables en GaAs n-type)
        phi = resolver_poisson(n, np.zeros_like(n), phi_0, phi_N, dx, param)

        # 4. criterio de convergencia
        delta = np.max(np.abs(n - n_prev)) / (np.max(n) + 1e-30)
        if delta < tol:
            conv = True
            break

    if not conv:
        print(f"  ⚠ No convergió en {max_iter} pasos (delta={delta:.2e})")

    return n, phi, calcular_corriente(n, phi, dx, param), conv


def barrido_tension(V_vals, n_init, nisc, n0_L, phi_N, dx, dt, param):
    """
    Resuelve el estado estacionario para cada V_app de V_vals, usando la
    solución del voltaje anterior como condición inicial (warm start).
    Imprime una tabla con los resultados.

    Retorna: J_vals [A/cm²] (array), n del último voltaje
    """
    n_actual = n_init
    J_vals   = []

    print(f"\n{'V_app':>8}  {'j [A/cm²]':>14}  {'conv':>5}")
    print("-" * 35)

    for V_app in V_vals:
        n_actual, _, j, conv = simular_vapp(
            V_app, n_actual, nisc, n0_L, phi_N, dx, dt, param
        )
        J_vals.append(j)
        print(f"{V_app:+8.3f}  {j:+14.4e}  {'✓' if conv else 'X'}")

    return np.array(J_vals), n_actual


def perfiles_en_tensiones(V_lista, n_base, nisc, n0_L, phi_N, dx, dt, param):
    """
    Perfiles estacionarios de n y phi para cada V_app de V_lista,
    todos partiendo de la misma condición inicial n_base.

    Retorna: lista de tuplas (n, phi), una por tensión
    """
    perfiles = []
    for V_app in V_lista:
        n, phi, _, _ = simular_vapp(V_app, n_base, nisc, n0_L, phi_N, dx, dt, param)
        perfiles.append((n, phi))
    return perfiles
