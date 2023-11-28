# Author: Loren Matilsky
# Created: 11/28/2023
#
# Purpose: generate a binary file (for Rayleigh to read) that contains
# a heating profile Q(r) = c_10*f_6
# that enforces a convectively unstable layer below a (slightly) 
# stable layer
#
# parameters:
# alpha: aspect ratio of WL-to-CZ
# beta: aspect ratio of full layer
# delta1: width of bottom heating layer
# delta2: width of top cooling layer
# fluxratio: ratio of stable flux to unstable flux

import numpy as np
import sys, os

sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])

from reference_tools import equation_coefficients
from common import *
from cla_util import *

nr_default = 5000
alpha_default = 0.25
beta_default = 0.9
delta1_default = 0.1
delta2_default = 0.1
fluxratio_default = 10.

def integrate_from_r0(integrand, rr, r0):
    # basic grid info
    nr = len(rr)
    ir0 = np.argmin(np.abs(rr - r0))
   
    # compute indefinite integral
    integral = np.zeros(nr)
    for ir in range(nr):
        integral[ir] = simps(integrand[ir:ir0 + 1], rr[ir:ir0 + 1])
        if ir <= ir0:
            integral[ir] *= -1
    return integral

def generate_heating_CZ_WL(nr=nr_default,  alpha=alpha_default, beta=beta_default, delta1=delta1_default, delta2=delta2_default, fluxratio=fluxratio_default):

    # compute radial locations and grid
    rin = beta/(1.-beta)
    rout = 1./(1.-beta)
    r0 = rout - alpha/(alpha+1.)
    rr = np.linspace(rout, rin, nr)

    # make heating shape for CZ
    ir0 = np.argmin(np.abs(rr-r0))
    heating_cz = np.exp(-(rr-rin)/delta1) - np.exp((r0-rr)/delta2)
    heating_cz[:ir0+1] = 0.

    # then heating shape for WL
    nr_wl = ir0 + 1
    heating_wl = np.zeros(nr)
    heating_wl[:ir0+1] = np.linspace(fluxratio, -fluxratio, nr_wl)

    heating = heating_cz + heating_wl
    return rin, rout, r0, rr, heating
