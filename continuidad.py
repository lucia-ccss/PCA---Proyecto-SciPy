# continuidad.py
import numpy as np
import parametros as param


def paso_continuidad(n, phi, G, dt, dx, v_e, n0_0, n0_L):
    """
    Avanza un paso temporal la ecuación de continuidad de electrones
    (esquema explícito: adelante en el tiempo, centrado en el espacio).

    Recibe:
        n    : array de largo N+1, densidad de electrones en el tiempo m [cm⁻³]
        phi  : array de largo N+1, potencial en el tiempo m [V]
        G    : array de largo N+1, tasa de generación [cm⁻³ s⁻¹]
        dt   : paso temporal [s]
        dx   : paso espacial [cm]
        v_e  : velocidad de recombinación superficial en Schottky [cm/s]
        n0_0 : densidad de equilibrio en x=0 [cm⁻³]
        n0_L : densidad de equilibrio en x=L [cm⁻³]

    Devuelve:
        n_nuevo : array de largo N+1, densidad de electrones en t+dt [cm⁻³]

    Usa de parametros.py: mu_e, D_e, tau_e

    Condiciones de contorno:
        x=0 (Schottky): -J_e(0) = v_e * (n_0 - n0_0)
        x=L (óhmico):   n_N = n0_L

    Estabilidad: requiere param.D_e * dt / dx**2 <= 1/2.
    """
    N = len(n) - 1

    # ── flujo J en los N semipuntos j+½  ──────────────────────────────────
    dphi  = (phi[1:] - phi[:-1]) / dx
    n_mid = 0.5 * (n[1:] + n[:-1])
    dn    = (n[1:]  - n[:-1])  / dx
    J     = param.mu_e * dphi * n_mid - param.D_e * dn

    # ── flujo entrante en cada nodo: J_izq[j] ≡ J_{j−½}  ─────────────────
    J_izq     = np.empty(N + 1)
    J_izq[0]  = -v_e * (n[0] - n0_0)      # BC Schottky
    J_izq[1:] = J

    # ── actualización explícita para j = 0 … N−1  ─────────────────────────
    G = np.broadcast_to(G, n.shape)
    n_nuevo       = np.empty(N + 1)
    n_nuevo[:-1]  = (n[:-1]
                     - dt / dx          * (J - J_izq[:-1])
                     - dt / param.tau_e * n[:-1]
                     + dt               * G[:-1])

    # ── BC óhmico en j = N  ───────────────────────────────────────────────
    n_nuevo[-1] = n0_L

    return n_nuevo