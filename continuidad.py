#continuidad

def paso_continuidad(n, phi, G, dt, dx, mu_e, D_e, tau_e, v_e, n0_0, n0_L):
    """
    Avanza un paso temporal la ecuación de continuidad de electrones
    (esquema explícito: adelante en el tiempo, centrado en el espacio).

    Recibe:
        n      : array de largo N+1, densidad de electrones en el tiempo m [cm⁻³]
        phi    : array de largo N+1, potencial en el tiempo m [V]
        G      : array de largo N+1, tasa de generación [cm⁻³ s⁻¹]
        dt     : paso temporal [s]
        dx     : paso espacial [cm]
        mu_e   : movilidad de electrones [cm²/(V·s)]
        D_e    : constante de difusión [cm²/s]
        tau_e  : tiempo de vida [s]
        v_e    : velocidad de recombinación superficial en el contacto Schottky [cm/s]
        n0_0   : densidad de equilibrio en x=0, n_e⁰(0) [cm⁻³]
        n0_L   : densidad de equilibrio en x=L, n_e⁰(L) [cm⁻³]

    Devuelve:
        n_nuevo : array de largo N+1, densidad de electrones en el tiempo m+1 [cm⁻³]

    Condiciones de contorno:
        x=0 (Schottky): -J_e(0) = v_e * (n_0 - n0_0)
        x=L (óhmico):   n_N = n0_L

    Estabilidad: requiere D_e*dt/dx**2 <= 1/2.
    """
    ...