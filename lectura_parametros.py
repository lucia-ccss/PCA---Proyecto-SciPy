# lectura_parametros.py
"""
Lectura de parámetros (parametros.txt) y constantes (constantes.txt).

Formato de los .txt:  nombre = valor   # comentario
El valor puede ser un número o una expresión aritmética simple
(ej: 0.063 * 9.11e-31). Las líneas vacías y las que empiezan con
'#' se ignoran.
"""
import ast
import operator
from types import SimpleNamespace

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


def cargar_parametros(*rutas):
    """
    Lee uno o más .txt y junta todo en un único objeto con los parámetros
    como atributos (param.mu_e, param.T, ...). Falla si un mismo nombre
    aparece en más de un archivo.
    """
    todos = {}
    for ruta in rutas:
        nuevos = leer_txt(ruta)
        repetidos = todos.keys() & nuevos.keys()
        if repetidos:
            raise ValueError(f"Parámetros definidos más de una vez: {sorted(repetidos)}")
        todos.update(nuevos)
    return SimpleNamespace(**todos)
