# Author: Loren Matilsky
# Created:t st
`

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
    return dict({'density': rho_c*rho_nd, 'pressure': P_c*P_nd,\
            'temperature': T_c*T_nd, 'entropy': S_c*np.ones(nr)})

def compute_polytrope2(Nrho=5, r0=rm, ro=None, rho0=rhom, T0=Tm, mstar=msun, poly_n=1.5, gas_constant_star=thermo_R, nr=5000):
    poly_a = G*mstar/((poly_n + 1)*gas_constant_star*T0*r0)
    if ro is None:
        rmax = r0*poly_a/(poly_a - 1) # at rmax, rho = 0
        rtmp = np.linspace(r0, rmax, 100000) # make this grid super fine
        rho_ratio = (poly_a*(r0/rtmp) + (1 - poly_a))**poly_n
        Nrho_tmp = np.log(1/rho_ratio)
        iro = np.argmin(np.abs(Nrho_tmp - Nrho))
        ro = rtmp[iro]

    r = np.linspace(ro, ri, nr)
    temp = T0*(poly_a*(r0/r) + (1 - poly_a))
    rho = rho0*(poly_a*(r0/r) + (1 - poly_a))**poly_n
    prs = rho*gas_constant_star*temp
    cv_star = gas_constant_star/(thermo_gamma - 1)
    dsdr = cv_star*(poly_n/1.5 - 1)/(r + (1 - poly_a)*r**2/(poly_a*r0))
    entropy = cv_star*(poly_n/1.5 - 1)*(np.log(r/r0) - np.log(poly_a + (1 - poly_a)*(r/r0)))
    return dict({'rho': rho, 'prs': prs, 'temp': temp, 'dsdr': dsdr, 'entropy': entropy})
