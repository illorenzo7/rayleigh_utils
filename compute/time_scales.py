# Author: Loren Matilsky
# Created: 01/31/2020
# This script contains routines for computing various time scales
# for the Rayleigh run in directory [dirname],

# compute_tdt(dirname, mag=False) computes the thermal
# (default; magnetic, if desired)
# diffusion time using the (constant) shell depth as the length scale. 
# Gets diffusion profiles from transport or equation_coefficients

# compute_Prot(dirname) computes the rotation period, getting 
# angular_velocity from either main_input or equation_coefficients

import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from get_parameter import get_parameter
from get_eq import get_eq

# Get directory name
dirname = sys.argv[1]

def compute_tdt(dirname, mag=False, visc=False, tach=False):
    # Returns computed diffusion time (in sec) across whole layer
    # If tach=True, return diffusion time across whole layer,
    # across CZ and across RZ (tuple of 3)
    # Read in the diffusion profile
    eq = get_eq(dirname)
    rr = eq.radius
    if mag:
        diff = eq.eta
    elif visc:
        diff = eq.nu
    else:
        diff = eq.kappa

    # Compute and return the diffusion time
    if tach:
        domain_bounds = get_parameter(dirname, 'domain_bounds')
        ri, rm, ro = domain_bounds
        rmid = 0.5*(ri + ro)
        rmidrz = 0.5*(ri + rm)
        rmidcz = 0.5*(rm + ro)

        irmidrz = np.argmin(np.abs(rr - rmidrz))
        irmidcz = np.argmin(np.abs(rr - rmidcz))
        irmid = np.argmin(np.abs(rr - rmid))

        diff_midrz = diff[irmidrz]
        diff_midcz = diff[irmidcz]
        diff_mid = diff[irmid]

        Hrz = rm - ri
        Hcz = ro - rm
        H = ro - ri

        return Hrz**2.0/diff_midrz, Hcz**2.0/diff_midcz, H**2.0/diff_mid
    else:
        ri, ro = np.min(rr), np.max(rr)
        rmid = 0.5*(ri + ro)
        irmid = np.argmin(np.abs(rr - rmid))
        diff_mid = diff[irmid]
        H = ro - ri
        return H**2.0/diff_mid

def compute_Prot(dirname):
    try:
        Om0 = get_parameter(dirname, 'angular_velocity')
    except:
        eq = equation_coefficients()
        eq.read(dirname + '/equation_coefficients')
        Om0 = eq.constants[0]/2.
    return 2*np.pi/Om0
