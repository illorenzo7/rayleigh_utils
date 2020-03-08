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
from rayleigh_diagnostics import TransportCoeffs
from reference_tools import equation_coefficients

# Get directory name
dirname = sys.argv[1]

def compute_tdt(dirname, mag=False, tach=False):
    # Returns computed diffusion time (in sec) across whole layer
    # If tach=True, return diffusion time across whole layer,
    # across CZ and across RZ (tuple of 3)
    # Read in the diffusion profile
    try: 
        trans = TransportCoeffs(dirname + '/transport')
        if mag:
            diff = trans.eta
        else:
            diff = trans.kappa
        rr = trans.radius
        print ("compute_tdt(): Got diffusion time from 'transport' file")
    except:
        eq = equation_coefficients()
        eq.read(dirname + '/equation_coefficients')
        if mag:
            diff = eq.constants[6]*eq.functions[6]
        else:
            diff = eq.constants[5]*eq.functions[4]
        rr = eq.radius
        print ("compute_tdt(): Got diffusion time from 'equation_coefficients' file")

    # Compute and return the diffusion time
    if tach:
        domain_bounds = get_parameter(dirname, 'domain_bounds')
        ri, rm, ro = domain_bounds
        iri = np.argmin(np.abs(rr - ri))
        irm = np.argmin(np.abs(rr - rm))
        iro = np.argmin(np.abs(rr - ro))
        diff_top_CZ = diff[iro]
        diff_top_RZ = diff[irm]
        H = ro - ri
        H_CZ = ro - rm
        H_RZ = rm - ri
        return H**2.0/diff_top_CZ, H_CZ**2.0/diff_top_CZ,\
                H_RZ**2.0/diff_top_RZ
    else:
        diff_top = diff[0]
        H = np.max(rr) - np.min(rr)
        return H**2.0/diff_top

def compute_Prot(dirname):
    try:
        Om0 = get_parameter(dirname, 'angular_velocity')
        print ("compute_Prot(): Got Prot from 'main_input' file")
    except:
        eq = equation_coefficients()
        eq.read(dirname + '/equation_coefficients')
        Om0 = eq.constants[0]/2.
        print ("compute_Prot(): Got Prot from 'equation_coefficients' file")
    return 2*np.pi/Om0
