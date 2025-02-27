# Author: Loren Matilsky
# Updated: 06/30/2021

import numpy as np
from scipy.integrate import simpson
import sys, os
sys.path.append(os.environ['raco'])
from common import g_univ, sun
from cla_util import *


kw_compute_polytrope2_default = dict({'nrho': 3, 'r0': sun.rbcz, 'r1': None, 'rho0': sun.rho_bcz, 'T0': sun.tmp_bcz, 'poly_mass': msun, 'poly_n': 1.5, 'gas_constant_star': sun.gas_constant})

def compute_polytrope2(**kw):
    # check for bad keys
    find_bad_keys(kw_compute_polytrope2_default, kw, 'compute_polytrope2()')
    # overwrite defaults
    kw = update_dict(kw_compute_polytrope2_default, kw)

    poly_a = g_univ*kw.poly_mass/((kw.poly_n + 1)*kw.gas_constant_star*kw.T0*kw.r0)
    nr = 10000 # make this grid super fine
    if kw.r1 is None:
        rmax = kw.r0*poly_a/(poly_a - 1) # at rmax, rho = 0
        rtmp = np.linspace(kw.r0, rmax, nr)
        rho_ratio = (poly_a*(kw.r0/rtmp) + (1 - poly_a))**kw.poly_n
        nrho_tmp = np.log(1/rho_ratio)
        ir1 = np.argmin(np.abs(nrho_tmp - kw.nrho))
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
    cv_star = kw.gas_constant_star/(gamma_ideal - 1)

    # S and its derivative will be zero most of the time...
    dSdr = cv_star*(kw.poly_n/1.5 - 1)/(r + (1 - poly_a)*r**2/(poly_a*kw.r0))
    S = cv_star*(kw.poly_n/1.5 - 1)*(np.log(r/kw.r0) - np.log(poly_a + (1 - poly_a)*(r/kw.r0)))

    return dict({'rr': r, 'rho': rho, 'dlnrho': dlnrho, 'd2lnrho': dlnrho, 'P': P, 'T': T,  'dlnT': dlnT, 'dSdr': dSdr, 'S': S})
