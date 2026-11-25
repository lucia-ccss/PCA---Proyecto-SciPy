# constantes del material

# parámetros de GaAs (según ioffe)
mu_e      = 8500                 # movilidad de los electrones [cm²/Vs]
D_e       = 200                  # coeficiente de difusión de los electrones [cm²/s]
tau_e     = 5e-9                 # tiempo de vida de los electrones [s]
e         = 1.602e-19            # carga del electrón [C]
epsilon_0 = 8.854e-14            # constante dieléctrica del vacío [F/cm]
epsilon_r = 12.9                 # constante dieléctrica relativa
m_eeff    = 0.063 * 9.11e-31     # masa efectiva e- [kg]
m_heff    = 0.51  * 9.11e-31     # masa efectiva h+ [kg]
E_g       = 1.424                # bandgap [eV]
C         = 1e14                 # dopaje [cm⁻³]

# otras constantes físicas
k_B  = 8.617333262e-5            # constante de Boltzmann [eV/K]
k_BJ = 1.380649e-23              # constante de Boltzmann [J/K]
h    = 6.626e-34                 # constante de Planck [J·s]

# parámetros del contacto GaAs–Au
phi_B = 0.7                      # altura de barrera Schottky (Fermi-pinning) [V]
v_e   = 1e7                      # velocidad de recombinación superficial [cm/s]

# parámetros condicionales
T = 300                          # temperatura [K]
L = 1e-4                         # longitud [cm] (un micrón)
N = 200                          # número de pasos