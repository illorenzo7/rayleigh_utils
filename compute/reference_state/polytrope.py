# Author: Loren Matilsky
# Updated: 06/30/2021

import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
from common import *
from cla_util import *

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
    return dict({'rho': rho_c*rho_nd, 'P': P_c*P_nd,\
            'T': T_c*T_nd, 'S': S_c*np.ones(nr)})

compute_polytrope2_kwargs_default = dict({'Nrho': 3, 'r0': rbcz, 'r1': None, 'rho0': rhobcz, 'T0': tempbcz, 'mstar': msun, 'poly_n': 1.5, 'gas_constant_star': thermo_R})

def compute_polytrope2(**kwargs):
    # overwrite defaults
    kw = update_dict(compute_polytrope2_kwargs_default, kwargs)
    # check for bad keys
    find_bad_keys(compute_polytrope2_kwargs_default, kwargs, 'compute_polytrope2()')

    poly_a = G*kw.mstar/((kw.poly_n + 1)*kw.gas_constant_star*kw.T0*kw.r0)
    nr = 100000 # make this grid super fine
    if kw.r1 is None:
        rmax = kw.r0*poly_a/(poly_a - 1) # at rmax, rho = 0
        rtmp = np.linspace(kw.r0, rmax, nr)
        rho_ratio = (poly_a*(kw.r0/rtmp) + (1 - poly_a))**kw.poly_n
        Nrho_tmp = np.log(1/rho_ratio)
        ir1 = np.argmin(np.abs(Nrho_tmp - kw.Nrho))
        kw.r1 = rtmp[ir1]

    r = np.linspace(kw.r1, kw.r0, nr) 
    # keep radii reversed in keeping with Rayleigh's convention

    factor = poly_a*kw.r0/r + (1 - poly_a)
    dfactor = -poly_a*kw.r0/r**2
    d2factor = 2*kw.poly_n*kw.r0/r**3
    T = kw.T0*factor
    dlnT = 1/factor*dfactor
    rho = kw.rho0*factor**kw.poly_n
    dlnrho = kw.poly_n/factor*dfactor
    d2lnrho = -kw.poly_n/factor**2 * dfactor**2*d2factor
    
    P = rho*kw.gas_constant_star*T
    cv_star = kw.gas_constant_star/(thermo_gamma - 1)

    # S and its derivative will be zero most of the time...
    dSdr = cv_star*(kw.poly_n/1.5 - 1)/(r + (1 - poly_a)*r**2/(poly_a*kw.r0))
    S = cv_star*(kw.poly_n/1.5 - 1)*(np.log(r/kw.r0) - np.log(poly_a + (1 - poly_a)*(r/kw.r0)))

    return dict({'rr': r, 'rho': rho, 'dlnrho': dlnrho, 'd2lnrho': dlnrho, 'P': P, 'T': T,  'dlnT': dlnT, 'dSdr': dSdr, 'S': S})
