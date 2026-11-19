# validacion_param.py
"""Chequeos de consistencia de parámetros antes de correr la simulación."""
import warnings
import numpy as np
from typing import Annotated
from pydantic import TypeAdapter, Field
from pydantic import ValidationError

# ── Adaptadores de tipo reutilizables ─────────────────────────────────────────
# TypeAdapter valida un valor suelto con las reglas de Pydantic,
# sin necesidad de definir una clase BaseModel.
_pos_float = TypeAdapter(Annotated[float, Field(gt=0)])   # float estrictamente positivo
_nn_float  = TypeAdapter(Annotated[float, Field(ge=0)])   # float no negativo
_int_ge3   = TypeAdapter(Annotated[int,   Field(ge=3)])   # entero >= 3


# ── Validación de parámetros físicos ─────────────────────────────────────────

def validar_fisicos(p):
    """Verifica que los atributos del módulo p tengan valores físicamente válidos."""
    campos_positivos = {
        "mu_e": p.mu_e, "D_e": p.D_e, "tau_e": p.tau_e,
        "epsilon_r": p.epsilon_r, "E_g": p.E_g,
        "m_eeff": p.m_eeff, "m_heff": p.m_heff,
        "T": p.T, "L": p.L, "k_B": p.k_B,
    }
    for nombre, valor in campos_positivos.items():
        try:
            _pos_float.validate_python(float(valor))
        except ValidationError:
            raise ValueError(f"{nombre} debe ser finito y > 0 (recibido {valor!r})")

    try:
        _nn_float.validate_python(float(p.C))
    except ValidationError:
        raise ValueError(f"C (dopaje) no puede ser negativo (recibido {p.C!r})")

    try:
        _int_ge3.validate_python(int(p.N))
    except ValidationError:
        raise ValueError(f"N debe ser un entero >= 3 (recibido {p.N!r})")


# ── Relación de Einstein ──────────────────────────────────────────────────────

def validar_einstein(p, rtol=0.1):
    """D_e debería ser ≈ mu_e · k_B · T. Avisa si la discrepancia supera rtol."""
    D_esperado = p.mu_e * p.k_B * p.T
    err = abs(p.D_e - D_esperado) / D_esperado
    if err > rtol:
        warnings.warn(
            f"D_e={p.D_e} difiere {err:.1%} de mu_e·k_B·T={D_esperado:.1f} cm²/s"
        )
    return err


# ── Estabilidad del esquema numérico ─────────────────────────────────────────

def validar_estabilidad(p, dt, mu=None, E_max=None):
    """
    Verifica D_e·dt/dx² ≤ 1/2.
    Opcionalmente revisa el número de Péclet de malla.
    Devuelve el valor de D_e·dt/dx².
    """
    _pos_float.validate_python(float(dt))   # dt también debe ser positivo

    dx = p.L / p.N
    r  = p.D_e * dt / dx**2
    if r > 0.5:
        raise ValueError(
            f"Esquema inestable: D_e·dt/dx² = {r:.3g} > 1/2. "
            f"Usá dt ≤ {0.5 * dx**2 / p.D_e:.3e} s"
        )
    if mu is not None and E_max is not None:
        Pe = mu * abs(E_max) * dx / (2 * p.D_e)
        if Pe > 1:
            warnings.warn(
                f"Péclet de malla = {Pe:.3g} > 1: "
                f"posibles oscilaciones; refiná dx o usá Scharfetter–Gummel"
            )
    return r


# ── Validación de arrays en tiempo de ejecución ───────────────────────────────

def validar_arrays(N, **arrays):
    """Uso: validar_arrays(N, n=n, phi=phi, G=G)"""
    for nombre, a in arrays.items():
        a = np.asarray(a)
        if a.shape != (N + 1,):
            raise ValueError(f"{nombre} tiene forma {a.shape}, se esperaba ({N+1},)")
        if not np.all(np.isfinite(a)):
            raise ValueError(f"{nombre} contiene NaN o inf")


# ── Punto de entrada único ────────────────────────────────────────────────────

def validar(p, dt):
    """Llama a todas las validaciones en orden. Uso: validar(param, dt)"""
    validar_fisicos(p)
    validar_einstein(p)
    return validar_estabilidad(p, dt)