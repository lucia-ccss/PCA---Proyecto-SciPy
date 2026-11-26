# Proyecto final

Modelar una juntura GaAs–Au como la siguiente:

```
      Au | GaAs | Au
         ↓      ↓
     contacto  contacto
       tipo      tipo
     Schottky   Óhmico
```

con el objetivo de calcular la corriente que se mide entre contactos para determinadas tensiones.

La descripción será unidimensional:

```
        x=0    GaAs    x=L
   Au  |////////////////|  Au
     Schottky         Óhmico
```

Queremos ver cómo se comportan los portadores de carga en el semiconductor.

## Ecuaciones de continuidad

$$
\partial_t n_e(x,t) = -\partial_x J_e(x,t) - R_e(x,t) + G_e(x,t)
$$

donde $n_e$ es la densidad de electrones, $J_e$ la densidad de corriente, $R_e$ la recombinación y $G_e$ la generación, con:

- $J_e(x,t) = \mu_e\,\partial_x \phi(x,t)\, n_e(x,t) - D_e\,\partial_x n_e(x,t)$

  donde $\mu_e$ es la movilidad electrónica y $D_e$ la constante de difusión de los portadores.

- $R_e(x,t) = \dfrac{n_e(x,t)}{\tau_e}$

  donde $\tau_e$ es el tiempo de vida.

La corriente estará dada por:

$$
j_e = -e J_e = -e\mu_e\,\partial_x \phi(x,t)\, n_e(x,t) + e D_e\,\partial_x n_e(x,t)
$$

y la densidad de carga por:

$$
\rho_e(x,t) = -e\, n_e(x,t)
$$

En resumen:

$$
\left[\begin{aligned}
\partial_t n_e &= -\partial_x J_e - R_e + G_e \\
J_e &= \mu_e\,\partial_x \phi\, n_e - D_e\,\partial_x n_e
\end{aligned}\right]
$$

$$
\Rightarrow \quad \partial_t n_e = -\mu_e\,\partial^2_{xx}\phi\; n_e - \mu_e\,\partial_x \phi\,\partial_x n_e + D_e\,\partial^2_{xx} n_e - R_e + G_e
$$

Y para huecos es todo análogo, con las constantes que correspondan.

## Ecuación de Poisson

La **ecuación de Poisson** que hay que resolver en el semiconductor es:

$$
-\partial_x\left[\varepsilon_r(x)\,\partial_x \phi(x)\right] = \frac{e}{\varepsilon_0}\left[n_h(x,t) - n_e(x,t) + C(x,t)\right]
$$

donde $C(x,t) = N_d^+(x,t) - N_a^-(x,t)$ es el perfil de impurezas.

## Condiciones de contorno

En todo el problema el potencial se relaciona con la densidad de electrones en equilibrio a través de la relación de Boltzmann

$$
n_e^0(x) = n_{isc}\,\exp\left[\frac{e\,\phi^0(x)}{k_B T}\right]
\quad\Longleftrightarrow\quad
\phi^0(x) = \frac{k_B T}{e}\,\ln\left[\frac{n_e^0(x)}{n_{isc}}\right],
$$

y esa misma convención se usa en los dos contactos.

**Contacto tipo Schottky ($x = 0$):**

$$
-J_e(0,t) = v_e\left[n_e(0,t) - n_e^0(0)\right]
$$

donde $v_e$ es la velocidad de recombinación superficial y $n_e^0(0)$ la densidad de portadores en el equilibrio. El signo menos es por convención. Una corriente $J_e$ positiva es el flujo hacia $+x$. Si $n_e(0) > n_e^0(0)$, los electrones salen hacia el metal y $J_e(0) < 0$. 

Por los estados de superficie del GaAs, el nivel de Fermi queda anclado (Fermi-pinning) en la interfaz. Eso fija la **altura de barrera** $\varphi_B = (E_c - E_F)/e \approx 0{,}7\ \mathrm{V}$, independientemente del metal. La densidad de electrones en equilibrio en la interfaz es entonces

$$
n_e^0(0) = N_c\,\exp\left[-\frac{e\,\varphi_B}{k_B T}\right],
\qquad
N_c = 2\left[\frac{2\pi m_e^{\ast}\, k_B T}{h^2}\right]^{3/2},
$$

con $N_c$ la densidad efectiva de estados de la banda de conducción. Usando la relación de Boltzmann, el potencial en el contacto es

$$
\phi(0,t) = \frac{k_B T}{e}\,\ln\left[\frac{n_e^0(0)}{n_{isc}}\right] + V_{App}
= \frac{k_B T}{e}\,\ln\left[\frac{N_c}{n_{isc}}\right] - \varphi_B + V_{App}.
$$

Para GaAs a 300 K esto da $n_e^0(0) \sim 10^{6}\ \mathrm{cm^{-3}}$ y $\phi(0) \approx -0{,}03\ \mathrm{V}$ en equilibrio, de modo que cerca del metal hay una zona de vaciamiento. $V_{App} > 0$ reduce la curvatura de bandas (polarización directa) y $V_{App} < 0$ la aumenta (polarización inversa).

**Contacto Óhmico ($x = L$):**

$$
n_e(L,t) = n_e^0(L)
$$

$$
\phi(L,t) = \varphi_{bi}(L)
$$

La tensión $V_{App}$ se aplica sólo en el contacto Schottky y el contacto óhmico queda como referencia. 

El potencial de referencia en el contacto óhmico es

$$
\varphi_{bi}(L) = \frac{k_B T}{e}\,\ln\left[\frac{n_e^0(L)}{n_{isc}}\right]
$$

con $n_e^0(L)$ la densidad de portadores en equilibrio:
$$
n_e^0(L) = \frac{C(L)}{2} + \sqrt{\frac{C(L)^2}{4} + n_{isc}^2}
$$

y $n_{isc}$ la densidad de portadores del semiconductor intrínseco:

$$
n_{isc} = 2\left[\frac{2\pi k_B T}{h^2}\right]^{3/2}\left(m_e^{\ast}\, m_h^{\ast}\right)^{3/4} e^{-E_g / 2 k_B T}.
$$

La diferencia $\varphi_{bi}(L) - \phi(0) \approx 0{,}67\ \mathrm{V}$ (para $C(L) = 10^{17}\ \mathrm{cm^{-3}}$ y $V_{App}=0$) es la barrera de potencial interna de la juntura.

## Diferencias finitas

Llevamos el problema a algo resoluble computacionalmente, planteando las ecuaciones en **diferencias finitas**.

Para eso, discretizamos el tiempo como $t = m\,\Delta t$ (con $m = 0, \dots, M$ y $\Delta t = t_{max}/M$) y el espacio como $x = j\,\Delta x$ (con $j = 0, \dots, N$ y $\Delta x = L/N$).

Notación: $n_j^m \equiv n_e(x_j, t_m)$, $p_j^m \equiv n_h(x_j, t_m)$, $\phi_j^m \equiv \phi(x_j, t_m)$, $C_j^m \equiv C(x_j, t_m)$ y $G_j^m \equiv G_e(x_j, t_m)$.

### Ecuación de continuidad

Usando diferencias hacia adelante en el tiempo y centradas en el espacio:

$$
\begin{aligned}
n_j^{m+1} = n_j^m
&- \mu_e \frac{\Delta t}{\Delta x^2}\, n_j^m \left[\phi_{j+1}^m - 2\phi_j^m + \phi_{j-1}^m\right] \\
&- \mu_e \frac{\Delta t}{4\,\Delta x^2} \left[\phi_{j+1}^m - \phi_{j-1}^m\right]\left[n_{j+1}^m - n_{j-1}^m\right] \\
&+ D_e \frac{\Delta t}{\Delta x^2}\left[n_{j+1}^m - 2 n_j^m + n_{j-1}^m\right] \\
&- \Delta t\, \frac{n_j^m}{\tau_e} + \Delta t\, G_j^m
\end{aligned}
$$

para $j = 1, \dots, N-1$. Al ser un esquema explícito, requiere $D_e\,\Delta t/\Delta x^2 \le 1/2$ para ser estable. Además, si la deriva domina ($\mu_e |\Delta\phi| \gtrsim 2 D_e$ entre nodos vecinos), el centrado en el espacio puede oscilar; en ese caso conviene discretizar el flujo con el esquema de Scharfetter–Gummel.

### Ecuación de Poisson

Con $\varepsilon_r$ constante y diferencias centradas, para cada nodo interior $j = 1, \dots, N-1$:

$$
-\phi_{j-1}^m + 2\phi_j^m - \phi_{j+1}^m = \Delta x^2\, \frac{e}{\varepsilon_0 \varepsilon_r}\left[p_j^m - n_j^m + C_j^m\right]
$$

Es un problema de contorno: $\phi_0^m$ y $\phi_N^m$ están fijados por los contactos, así que no se puede avanzar nodo a nodo desde un extremo (no se impondría la condición del otro extremo y los errores se amplifican). Pasando los valores de borde al lado derecho queda un sistema tridiagonal de $N-1$ incógnitas:

$$
\begin{pmatrix}
2 & -1 & & \\
-1 & 2 & \ddots & \\
& \ddots & \ddots & -1 \\
& & -1 & 2
\end{pmatrix}
\begin{pmatrix} \phi_1^m \\ \phi_2^m \\ \vdots \\ \phi_{N-1}^m \end{pmatrix}
= \Delta x^2\, \frac{e}{\varepsilon_0 \varepsilon_r}
\begin{pmatrix} p_1^m - n_1^m + C_1^m \\ p_2^m - n_2^m + C_2^m \\ \vdots \\ p_{N-1}^m - n_{N-1}^m + C_{N-1}^m \end{pmatrix}
+ \begin{pmatrix} \phi_0^m \\ 0 \\ \vdots \\ \phi_N^m \end{pmatrix}
$$

que se resuelve en cada paso de tiempo con `scipy.linalg.solve_banded`.

### Condiciones de contorno

$$
\begin{aligned}
n_e(0,t) &\longrightarrow n_0^m \\
n_e(L,t) &\longrightarrow n_N^m \\
\phi(0,t) &\longrightarrow \phi_0^m \\
\phi(L,t) &\longrightarrow \phi_N^m
\end{aligned}
$$

#### Potencial 

Es una condición tipo Dirichlet. Los dos valores de borde son fijos para cada tensión aplicada.

$$
\phi_0^m = \frac{k_B T}{e}\,\ln\left[\frac{N_c}{n_{isc}}\right] - \varphi_B + V_{App},
\qquad \varphi_B \approx 0{,}7\ \mathrm{V}
$$

$$
\phi_N^m = \varphi_{bi}(L) = \frac{k_B T}{e}\,\ln\left[\frac{n_e^0(L)}{n_{isc}}\right]
$$

Si $V_{App}$ es constante durante la simulación, $\phi_0^m$ y $\phi_N^m$ no dependen de $m$.

#### Densidad en el contacto óhmico 
Condición tipo Dirichlet.

$$
n_N^m = n_e^0(L) = \frac{C(L)}{2} + \sqrt{\frac{C(L)^2}{4} + n_{isc}^2}
$$

#### Densidad en el contacto Schottky

Condición mixta Dirichlet y Neumann. Aproximando las derivadas de $J_e(0,t)$ con diferencias hacia adelante,

$$
J_e(0) = \mu_e\,\partial_x \phi(x,t)\, n_e(x,t) - D_e\,\partial_x n_e(x,t) \approx \mu_e\,\frac{\phi_1 - \phi_0}{\Delta x}\, n_0 - D_e\,\frac{n_1 - n_0}{\Delta x},
$$

y reemplazando en $-J_e(0,t) = v_e\left[n_e(0,t) - n_e^0(0)\right]$, se despeja $n_0$:

$$
n_0^{m+1} = \frac{v_e\, n_e^0(0) + \dfrac{D_e}{\Delta x}\, n_1^{m+1}}{v_e + \dfrac{D_e}{\Delta x} + \dfrac{\mu_e}{\Delta x}\left(\phi_1^m - \phi_0^m\right)}
$$

Se evalúa después de actualizar los nodos interiores, usando el $n_1^{m+1}$ ya calculado. La densidad de equilibrio en el contacto Schottky está fijada por la altura de barrera:

$$
n_e^0(0) = N_c\, \exp\left[-\frac{e\,\varphi_B}{k_B T}\right],
\qquad
N_c = 2\left[\frac{2\pi m_e^{\ast}\, k_B T}{h^2}\right]^{3/2}
$$

### Condición inicial

Se parte de la densidad de equilibrio $n_e^0(x)$ (sin generación ni tensión aplicada):

$$
n_e(x,0) = n_e^0(x) \;\longrightarrow\; n_j^0 = n_e^0(x_j)
$$
