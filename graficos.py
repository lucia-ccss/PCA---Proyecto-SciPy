# graficos.py
"""Gráficos de la curva I-V y de los perfiles de portadores / potencial."""
import matplotlib.pyplot as plt


def graficar_curva_IV(V_vals, J_vals, archivo="curva_IV.png"):
    """Curva I-V lineal, guardada en `archivo`."""
    fig, ax = plt.subplots()
    ax.plot(V_vals, J_vals * 1e3, 'o-', lw=1.5, ms=4)
    ax.axhline(0, color='k', lw=0.5)
    ax.axvline(0, color='k', lw=0.5)
    ax.set_xlabel("$V_{App}$ [V]")
    ax.set_ylabel("$j_e$ [mA/cm²]")
    ax.set_title("Curva I–V: juntura GaAs–Au (Schottky–Óhmico)")
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    fig.savefig(archivo, dpi=150)
    print(f"\nGuardada: {archivo}")
    return fig


def graficar_perfiles_y_IV(x_um, V_perfiles, perfiles, V_vals, J_vals,
                           archivo="perfiles_y_IV.png"):
    """
    Figura 2x2:
      - izquierda: perfiles de n_e (log) y phi para cada tensión de V_perfiles
      - derecha:   curva I-V en escala lineal y logarítmica

    perfiles: lista de tuplas (n, phi), una por cada valor de V_perfiles.
    """
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=False)

    # ── perfiles ─────────────────────────────────────────────────────
    for V_app, (n, phi) in zip(V_perfiles, perfiles):
        label = f"$V_{{App}}={V_app:g}$ V"
        axes[0, 0].semilogy(x_um, n, label=label)
        axes[1, 0].plot(x_um, phi, label=label)

    axes[0, 0].set_ylabel("$n_e$ [cm⁻³]")
    axes[0, 0].set_title("Densidad de electrones")
    axes[1, 0].set_ylabel("$\\phi$ [V]")
    axes[1, 0].set_title("Potencial electrostático")
    for ax in axes[:, 0]:
        ax.set_xlabel("$x$ [µm]")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.4)

    # ── curva I-V lineal ─────────────────────────────────────────────
    axes[0, 1].plot(V_vals, J_vals * 1e3, 'o-', lw=1.5, ms=4, color='C1')
    axes[0, 1].axhline(0, color='k', lw=0.5)
    axes[0, 1].set_ylabel("$j_e$ [mA/cm²]")
    axes[0, 1].set_title("Curva I–V (lineal)")
    axes[0, 1].set_xlabel("$V_{App}$ [V]")
    axes[0, 1].grid(True, alpha=0.4)

    # ── curva I-V log (|j| separando signos) ─────────────────────────
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
    fig.savefig(archivo, dpi=150)
    print(f"Guardada: {archivo}")
    return fig


def mostrar():
    """Muestra todas las figuras abiertas."""
    plt.show()
