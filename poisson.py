#poisson

def resolver_poisson(n, p, C, phi0, phiN, dx, eps_r):
    """
    Resuelve Poisson en los nodos interiores.

    Recibe:
        n, p   : arrays de largo N+1, densidades de electrones y huecos [cm⁻³]
        C      : array de largo N+1, perfil de dopaje [cm⁻³]
        phi0   : potencial en x=0 [V]
        phiN   : potencial en x=L [V]
        dx     : paso espacial [cm]
        eps_r  : permitividad relativa

    Devuelve:
        phi    : array de largo N+1 con el potencial en toda la grilla [V]
    """
    ...  # esto lo implementa quien hace Poisson