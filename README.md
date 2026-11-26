# PCA---Proyecto-SciPy
# PHOTOCONDUCTIVE ANTENNA
✨ Autoras ✨

Lucía Cervantes Schamun

Candela Lopez Porto

## Simulación de juntura GaAs–Au (Schottky–Óhmico)

Este proyecto evalúa el modelo unidimensional de un semiconductor III-V (GaAs) con dos contactos metálicos de Au: un contacto tipo Schottky en $x = 0$ y un contacto Óhmico en $x = L$.
El objetivo es calcular la corriente entre contactos en función de la tensión aplicada,
siguiendo la dinámica de portadores de carga (electrones y huecos) en el semiconductor.

### Modelo físico

- **Ecuaciones de continuidad** para electrones y huecos, con recombinación tipo Shockley–Read–Hall ($\tau_e$) y término de generación $G$.
- **Ecuación de Poisson** acoplada, con perfil de impurezas $C(x) = N_d^+ - N_a^-$.
- **Condición de contorno Schottky** ($x = 0$): velocidad de recombinación superficial $v_e$
  y potencial fijado por Fermi-Pinning ($\varphi_\text{Stky} = 0{,}7\ \text{V}$) más la tensión aplicada $V_\text{App}$.
- **Condición de contorno Óhmico** ($x = L$): densidad fijada en equilibrio $n_e^0(L)$
  y potencial built-in $\varphi_{bi}(L)$ como referencia.

### Método numérico

Las ecuaciones se discretizan en una grilla espacio-temporal mediante diferencias finitas:

- Esquema explícito (adelante en tiempo, centrado en espacio) para la ecuación de continuidad,   sujeto a la condición de estabilidad $D_e\,\Delta t/\Delta x^2 \le 1/2$.
  Para regímenes dominados por deriva se recomienda el esquema de Scharfetter–Gummel.
- Sistema tridiagonal para la ecuación de Poisson en cada paso de tiempo.

La condición inicial es la densidad de equilibrio $n_e^0(x)$ sin tensión aplicada.

Para ver el detalle de las ecuaciones utilizados y modelos considerados ingresar a [`PLANTEO.md`](PLANTEO.md).

Para ver cómo está organizado el código y qué hace cada archivo, ver [`README2.md`](README2.md).

Para este proyecto se deben instalar: numpy scipy matplotlib pydantic
o utilizar uv