# Estructura del código

## Qué hace el proyecto

Simula en 1D una juntura de GaAs con dos contactos de oro: un contacto
Schottky en $x=0$ y uno óhmico en $x=L$. Para cada tensión aplicada $V_{App}$
resuelve juntas la ecuación de continuidad de electrones y la de Poisson,
avanzando en el tiempo hasta llegar al estado estacionario. Con eso calcula
la corriente y arma la curva I-V.

El desarrollo de las ecuaciones está en [PLANTEO.md](PLANTEO.md).

## Cómo correrlo

```
uv run python main.py                  # modo interactivo: muestra T, L, N y pregunta si cambiarlos
uv run python main.py --T 350          # modo CLI: cambia T sin preguntar
uv run python main.py --T 350 --L 2e-4 --N 300
```

La primera vez, `uv sync` crea el entorno `.venv` con las dependencias.

## Cómo está construido `main.py`

`main.py` solo organiza el flujo del programa. Todas las funciones están en
los otros módulos. Los pasos de `main()` son:

1. **Parámetros**: `cargar_parametros("constantes.txt", "parametros.txt")`
   devuelve el objeto `param` (`param.mu_e`, `param.T`, ...). Ese objeto se
   pasa como argumento a todas las funciones.
2. **Configuración**: `configurar_parametros(param)` permite cambiar T, L y N
   por línea de comandos o de forma interactiva.
3. **Grilla**: `dx = L/N` y `dt = 0.4·dx²/D_e`, que deja $D\,\Delta t/\Delta x^2 = 0.4$,
   por debajo del límite de estabilidad de 0.5.
4. **Validación**: `validar(param, dt)` revisa que los parámetros sean
   físicos, la relación de Einstein y la estabilidad del esquema.
5. **Magnitudes de equilibrio**: `nisc` (densidad intrínseca), `n0_L` (densidad
   de equilibrio en el contacto óhmico) y `phi_N` (potencial del contacto
   óhmico).
6. **Equilibrio inicial**: `simular_vapp` a V=0, partiendo de una densidad
   uniforme `n0_L`.
7. **Barrido de tensión**: `barrido_tension` resuelve el estado estacionario
   de V = −1 a 0.5 V en 31 puntos. Cada tensión arranca desde la solución
   de la anterior (*warm start*), así converge más rápido.
8. **Perfiles**: `perfiles_en_tensiones` calcula $n(x)$ y $\phi(x)$ en
   V = 0, +0.3 y −0.5 V.
9. **Gráficos**: guarda `curva_IV.png` y `perfiles_y_IV.png` y muestra las
   figuras.

La función central es `simular_vapp` (en `simulacion.py`). Para una tensión
dada repite este ciclo:

1. avanza un paso de continuidad (esquema explícito),
2. corrige `n[0]` con la condición de contorno del Schottky (tipo Robin),
3. resuelve Poisson con p=0 (huecos despreciables en GaAs tipo n),
4. corta cuando el cambio relativo de `n` entre pasos es menor que `tol`.

Al final devuelve `n`, `phi`, la corriente y si convergió.

## Para qué sirve cada archivo

| Archivo | Para qué sirve |
|---|---|
| `main.py` | Programa principal (ver arriba). |
| `parametros.txt` | Parámetros del material GaAs: `mu_e`, `D_e`, `tau_e`, ε, masas efectivas, `E_g`, dopaje `C`, `phi_B`, `v_e`. |
| `parametrosSi.txt` | Los mismos parámetros para Si. |
| `constantes.txt` | Constantes físicas (`k_B`, `h`) y parámetros de la simulación (`T`, `L`, `N`). |
| `lectura_parametros.py` | Lee los `.txt` con formato `nombre = valor  # comentario`. Acepta expresiones como `0.063 * 9.11e-31` y las evalúa de forma segura con `ast`, sin usar `eval`. `cargar_parametros(...)` devuelve el objeto `param` y da error si un nombre aparece en más de un archivo. |
| `configurar_parametros.py` | Permite cambiar T, L y N con `argparse` (`--T`, `--L`, `--N`) o, si no se pasan argumentos, de forma interactiva, reintentando hasta recibir valores válidos. |
| `validacion_param.py` | Revisa los parámetros antes de correr, usando Pydantic: que sean positivos, que se cumpla la relación de Einstein $D \approx \mu k_B T$ (avisa si no) y que el esquema sea estable ($D\,\Delta t/\Delta x^2 \le 1/2$). También tiene `validar_arrays` para revisar forma y NaN de los arrays. |
| `equilibrio.py` | `calcular_nisc`, `calcular_Nc` y `n_equilibrio`. |
| `condiciones_contorno.py` | Potenciales en los contactos óhmico y Schottky, `n_eq_schottky` y `bc_n_schottky` (condición tipo Robin para `n[0]`). |
| `continuidad.py` | `paso_continuidad`: un paso temporal explícito de la ecuación de continuidad, con el flujo calculado en los puntos $j+1/2$. |
| `poisson.py` | `resolver_poisson`: arma el sistema tridiagonal y lo resuelve con `scipy.linalg.solve_banded`. |
| `simulacion.py` | `calcular_corriente`, `simular_vapp`, `barrido_tension` (el loop de la I-V) y `perfiles_en_tensiones`. |
| `graficos.py` | `graficar_curva_IV`, `graficar_perfiles_y_IV` y `mostrar`. |
| `PLANTEO.md` | Desarrollo de las ecuaciones y del modelo. |
| `bibliografía/` | PDFs de referencia sobre junturas y antenas fotoconductivas (PCA). |
| `pyproject.toml`, `uv.lock` | Dependencias del proyecto para `uv` (`uv sync` arma el `.venv`). |
| `src/pca_proyecto_scipy/` | Plantilla que generó `uv init`. No se usa todavía. |
| `curva_IV.png`, `perfiles_y_IV.png` | Gráficos que genera `main.py`. |
| `.gitignore` | Hace que git ignore `__pycache__/` y `.venv/`. |

## Dependencias

numpy, scipy, matplotlib, pydantic
