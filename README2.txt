ESTRUCTURA DEL CÓDIGO
=====================

Qué hace el proyecto
--------------------
Simula en 1D una juntura de GaAs con dos contactos de oro: un contacto
Schottky en x=0 y uno óhmico en x=L. Para cada tensión aplicada V_app
resuelve juntas la ecuación de continuidad de electrones y la de Poisson,
avanzando en el tiempo hasta llegar al estado estacionario. Con eso
calcula la corriente y arma la curva I-V.


Cómo está construido main.py
----------------------------
main.py solo organiza el flujo del programa. Todas las funciones están en
los otros módulos. Los pasos de main() son:

  1. Parámetros: cargar_parametros("constantes.txt", "parametros.txt")
     devuelve el objeto `param` (param.mu_e, param.T, ...). Ese objeto se
     pasa como argumento a todas las funciones.

  2. Grilla: dx = L/N y dt = 0.4·dx²/D_e, que deja D·dt/dx² = 0.4, por
     debajo del límite de estabilidad de 0.5.

  3. Validación: validar(param, dt) revisa que los parámetros sean
     físicos, la relación de Einstein y la estabilidad del esquema.

  4. Magnitudes de equilibrio: nisc (densidad intrínseca), n0_L (densidad
     de equilibrio en el contacto óhmico) y phi_N (potencial del contacto
     óhmico).

  5. Equilibrio inicial: simular_vapp a V=0, partiendo de una densidad
     uniforme n0_L.

  6. Barrido de tensión: barrido_tension resuelve el estado estacionario
     de V = -1 a 0.5 V en 31 puntos. Cada tensión arranca desde la
     solución de la anterior (warm start), así converge más rápido.

  7. Perfiles: perfiles_en_tensiones calcula n(x) y phi(x) en
     V = 0, +0.3 y -0.5 V.

  8. Gráficos: guarda curva_IV.png y perfiles_y_IV.png y muestra las
     figuras.

La función central es simular_vapp (en simulacion.py). Para una tensión
dada repite este ciclo:
  a. avanza un paso de continuidad (esquema explícito),
  b. corrige n[0] con la condición de contorno del Schottky (tipo Robin),
  c. resuelve Poisson con p=0 (huecos despreciables en GaAs tipo n),
  d. corta cuando el cambio relativo de n entre pasos es menor que tol.
Al final devuelve n, phi, la corriente y si convergió.


Para qué sirve cada archivo
---------------------------
parametros.txt
    Parámetros del material GaAs: mu_e, D_e, tau_e, epsilon, masas
    efectivas, E_g, dopaje C, phi_B, v_e.

constantes.txt
    Constantes físicas (k_B, h) y parámetros de la simulación (T, L, N).

lectura_parametros.py
    Lee los .txt con formato "nombre = valor  # comentario". Acepta
    expresiones como 0.063 * 9.11e-31 y las evalúa de forma segura con
    ast, sin usar eval. cargar_parametros(...) devuelve el objeto param y
    da error si un nombre aparece en más de un archivo.

validacion_param.py
    Revisa los parámetros antes de correr, usando Pydantic: que sean
    positivos, que se cumpla la relación de Einstein D ≈ mu·kT (avisa si
    no) y que el esquema sea estable (D·dt/dx² ≤ 1/2). También tiene
    validar_arrays para revisar forma y NaN de los arrays.

equilibrio.py
    calcular_nisc, calcular_Nc y n_equilibrio.

condiciones_contorno.py
    Potenciales en los contactos óhmico y Schottky, n_eq_schottky y
    bc_n_schottky (condición tipo Robin para n[0]).

continuidad.py
    paso_continuidad: un paso temporal explícito de la ecuación de
    continuidad, con el flujo calculado en los puntos j+1/2.

poisson.py
    resolver_poisson: arma el sistema tridiagonal y lo resuelve con
    scipy.linalg.solve_banded.

simulacion.py
    calcular_corriente, simular_vapp, barrido_tension (el loop de la I-V)
    y perfiles_en_tensiones.

graficos.py
    graficar_curva_IV, graficar_perfiles_y_IV y mostrar.

main.py
    Programa principal (ver arriba). Se corre con:  python main.py

proyecto_final.ipynb
    Notebook con el desarrollo de las ecuaciones y del modelo.

bibliografía/
    PDFs de referencia sobre junturas y antenas fotoconductivas (PCA).

pyproject.toml, src/pca_proyecto_scipy/
    Plantilla que generó "uv init". No se usa todavía y la lista de
    dependencias está vacía.

curva_IV.png, perfiles_y_IV.png
    Gráficos que genera main.py.

__pycache__/
    Archivos compilados de Python. No deberían estar en git; conviene
    agregarlos a un .gitignore.


Dependencias
------------
numpy, scipy, matplotlib, pydantic
