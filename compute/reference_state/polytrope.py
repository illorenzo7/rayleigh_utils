# Author: Loren Matilsky
# Created: well before 05/06/2019

# Computes the polytropic thermodynamic stratification for an 
# adiabatic atmosphere in a portion [ri, ro] of the solar CZ
# after the formulation of Jones et al. (2011).

import numpy as np
import sys, os
sys.path.append(os.envrione['raco'])
from common import *

def compute_polytrope(ri, ro, Nrho, nr, poly_n, rho_i):
    d = ro - ri
    beta = ri/ro
    poly_gamma = (poly_n + 1.)/poly_n
    msun = 1.98891e33
    gas_constant = c_P*(1. - 1./poly_gamma)

    r = np.linspace(ro, ri, nr)

    exp = np.exp(Nrho/poly_n)

    c0 = (1.+beta)/(1.-beta) * (1-beta*exp)/(1.+beta*exp)
    c1 = (1.+beta)*beta/(1.-beta)**2  * (exp - 1.)/(beta*exp + 1.)

    zeta = c0 + c1*d/r
    zeta_i = (1. + beta)*exp/(1. + beta*exp)

    rho_c = rho_i/zeta_i**poly_n

    T_c = G*msun/(c_P*c1*d)

    P_c = gas_constant*rho_c*T_c

    S_c = c_P*np.log(P_c**(1./poly_gamma)/rho_c) 
    # for an adiabatic polytrope,
            # S is constant everywhere
    rho_nd = zeta**(poly_n)
    P_nd = zeta**(poly_n+1.)
    T_nd = zeta
    return {'density': rho_c*rho_nd, 'pressure': P_c*P_nd,\
            'temperature': T_c*T_nd, 'entropy': S_c*np.ones(nr)}
