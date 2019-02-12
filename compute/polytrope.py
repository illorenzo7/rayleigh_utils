import numpy as np
import sys

# Compute the polytropic thermodynamic stratification for an 
# adiabatic atmosphere in a portion [ri, ro] of the solar CZ
def compute_polytrope(ri, ro, Nrho, nr, poly_n, rho_i):
    d = ro - ri
    beta = ri/ro
    poly_gamma = (poly_n + 1.)/poly_n
    cP = 3.5e8
    Newton_G = 6.67408e-8
    M_sun = 1.98891e33
    gas_constant = cP*(1. - 1./poly_gamma)

    r = np.linspace(ro, ri, nr)

    exp = np.exp(Nrho/poly_n)

    c0 = (1.+beta)/(1.-beta) * (1-beta*exp)/(1.+beta*exp)
    c1 = (1.+beta)*beta/(1.-beta)**2  * (exp - 1.)/(beta*exp + 1.)

    zeta = c0 + c1*d/r
    zeta_i = (1. + beta)*exp/(1. + beta*exp)

    rho_c = rho_i/zeta_i**poly_n

    T_c = Newton_G*M_sun/(cP*c1*d)

    P_c = gas_constant*rho_c*T_c

    S_c = cP*np.log(P_c**(1./poly_gamma)/rho_c) 
    # for an adiabatic polytrope,
            # S is constant everywhere
    rho_nd = zeta**(poly_n)
    P_nd = zeta**(poly_n+1.)
    T_nd = zeta
    return {'density': rho_c*rho_nd, 'pressure': P_c*P_nd,\
            'temperature': T_c*T_nd, 'entropy': S_c*np.ones(nr)}
