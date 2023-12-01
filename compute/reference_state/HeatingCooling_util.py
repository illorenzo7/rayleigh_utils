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
from scipy.integrate import simps
import sys, os

sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])

from reference_tools import equation_coefficients
from common import *
from cla_util import *

nr_default = 5000
alpha_default = 0.25
beta_default = 0.9
deltain_default = 0.1
deltaout_default = 0.1
deltac_default = 0.1
deltah_default = 0.1
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

# define quartic functions, each use heating layers
def psi_plus(rr, rc, delta):
    nr = len(rr)
    shape = np.zeros(nr)
    for ir in range(nr):
        rloc = rr[ir]
        if rloc <= rc:
            shape[ir] = 0.
        elif rloc < rc + delta:
            shape[ir] = (1. - ((rloc-rc)/delta)**2)**2
        else:
            shape[ir] = 0.
    return shape

def psi_minus(rr, rc, delta):
    return psi_plus(-(rr-2*rc), rc, delta)

def compute_heating_CZ_only(nr=nr_default, beta=beta_default, deltah=deltah_default, deltac=deltac_default):

    # compute radial locations and grid
    rin = beta/(1.-beta)
    rout = 1./(1.-beta)
    rr = np.linspace(rout, rin, nr)

    shape1 = psi_plus(rr, rin, deltah)
    shape2 = psi_minus(rr, rout, deltac)

    fourpi = 4*np.pi
    A1 = -1. / simps(fourpi*rr**2*shape1, rr) # remember rr is in decreasing order
    A2 = -1. / simps(fourpi*rr**2*shape2, rr)

    heating = A1*shape1 - A2*shape2

    return rin, rout, rr, heating

def compute_heating_CZ_WL(nr=nr_default, alpha=alpha_default, beta=beta_default, deltah1=deltah_default, deltac=deltac_default, deltah2=deltah_default, fluxratio=fluxratio_default):

    # compute radial locations and grid
    rin = beta/(1.-beta)
    r0 = 1./(1.-beta)
    rout = r0 + alpha
    rr = np.linspace(rout, rin, nr)

    shape1 = psi_plus(rr, rin, deltah1)
    shape2 = psi_minus(rr, r0, deltac) + psi_plus(rr, r0, deltac)
    shape3 = psi_minus(rr, rout, deltah2)

    fourpi = 4*np.pi
    A1 = -1. / simps(fourpi*rr**2*shape1, rr) # remember rr is in decreasing order
    A2 = -(1. + fluxratio) / simps(fourpi*rr**2*shape2, rr)
    A3 = -fluxratio / simps(fourpi*rr**2*shape3, rr)

    heating = A1*shape1 - A2*shape2 + A3*shape3

    return rin, r0, rout, rr, heating
