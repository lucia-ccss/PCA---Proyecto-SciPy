# configurar_parametros.py
"""
Gestión interactiva y por CLI de los parámetros condicionales de la simulación.
Uso:
    python main.py              → modo interactivo (pregunta con retry)
    python main.py --T 350      → modo CLI (aplica directo, falla si inválido)
    python main.py --T 350 --L 2e-4 --N 300
"""

import argparse


# ════════════════════════════════════════════════════════════════════
# Funciones de validación para argparse
# ════════════════════════════════════════════════════════════════════

def _tipo_T(valor):
    """Valida temperatura: float > 0."""
    try:
        v = float(valor)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"T debe ser un número (recibido '{valor}')"
        )
    if v <= 0:
        raise argparse.ArgumentTypeError(
            f"T debe ser > 0 (recibido {v})"
        )
    return v


def _tipo_L(valor):
    """Valida longitud: float > 0."""
    try:
        v = float(valor)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"L debe ser un número (recibido '{valor}')"
        )
    if v <= 0:
        raise argparse.ArgumentTypeError(
            f"L debe ser > 0 (recibido {v})"
        )
    return v


def _tipo_N(valor):
    """Valida número de nodos: entero >= 3."""
    try:
        v = int(valor)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"N debe ser un entero (recibido '{valor}')"
        )
    if v < 3:
        raise argparse.ArgumentTypeError(
            f"N debe ser >= 3 (recibido {v})"
        )
    return v


# ════════════════════════════════════════════════════════════════════
# Retry interactivo
# ════════════════════════════════════════════════════════════════════

def _pedir_con_retry(nombre, valor_actual, tipo_fn):
    """
    Solicita un valor por teclado y reintenta hasta recibir uno válido.
    Presionar Enter sin escribir nada conserva el valor actual.
    """
    while True:
        raw = input(f"  {nombre} [{valor_actual}]: ").strip()
        if not raw:
            return valor_actual
        try:
            return tipo_fn(raw)
        except argparse.ArgumentTypeError as e:
            print(f"  ✗ {e}. Intentá de nuevo.")


# ════════════════════════════════════════════════════════════════════
# Función principal del módulo
# ════════════════════════════════════════════════════════════════════

def configurar_parametros(param):
    """
    Configura T, L y N antes de la simulación.

    Recibe el objeto `param` leído de los .txt (cargar_parametros).

    · Si se pasan argumentos por CLI (--T, --L, --N), los aplica
      directamente con validación de argparse.
    · Si no se pasa ninguno, entra en modo interactivo con retry
      hasta que los valores sean válidos.

    Modifica param.T, param.L y param.N en el lugar.
    """
    parser = argparse.ArgumentParser(
        description="Simulación juntura GaAs–Au (Schottky–Óhmico)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--T", type=_tipo_T,
                        help="Temperatura [K]",
                        default=None)
    parser.add_argument("--L", type=_tipo_L,
                        help="Longitud del dispositivo [cm]",
                        default=None)
    parser.add_argument("--N", type=_tipo_N,
                        help="Número de nodos espaciales (>= 3)",
                        default=None)

    args = parser.parse_args()
    modo_cli = any(v is not None for v in vars(args).values())

    # ── Modo CLI ─────────────────────────────────────────────────────
    if modo_cli:
        if args.T is not None: param.T = args.T
        if args.L is not None: param.L = args.L
        if args.N is not None: param.N = args.N

        print("=" * 45)
        print("  Parámetros recibidos por CLI")
        print("=" * 45)
        print(f"  T = {param.T}    [K]")
        print(f"  L = {param.L}  [cm]")
        print(f"  N = {param.N}      [nodos]")
        print("=" * 45 + "\n")

    # ── Modo interactivo ─────────────────────────────────────────────
    else:
        print("=" * 45)
        print("  Parámetros de simulación")
        print("=" * 45)
        print(f"  T = {param.T}    [K]")
        print(f"  L = {param.L}  [cm]")
        print(f"  N = {param.N}      [nodos]")
        print("=" * 45)

        respuesta = input(
            "¿Continuar? [Enter para confirmar / 'no' para modificar]: "
        ).strip()

        if respuesta.lower() == "no":
            print("\nIngresá los nuevos valores (Enter para conservar el actual):\n")
            param.T = _pedir_con_retry("T", param.T, _tipo_T)
            param.L = _pedir_con_retry("L", param.L, _tipo_L)
            param.N = _pedir_con_retry("N", param.N, _tipo_N)
            print(
                f"\n  Valores actualizados: "
                f"T={param.T} K  |  L={param.L} cm  |  N={param.N}"
            )

        print("=" * 45 + "\n")

    return param