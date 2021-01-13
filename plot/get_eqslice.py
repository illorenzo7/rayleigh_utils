# Author: Loren Matilsky
# Created: 04/24/2020
#
# Extremely long and boring script to find fundamental fluid quantities
# or derivative quantities from an equatorial slice. 
# Takes an equatorial slice [eq] and
# varname in [vr, vt, vp, vl, vz, ...]
# returns the slice for the variable as an array of shape 
# (nphi, nr) 
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from varprops import var_indices, var_indices_old
from common import *

def get_eqslice(eq, varname, dirname=None, az=None, old=False, j=0):
    # Given an Equatorial_Slices object (nphi, ntheta, nr, nq, niter), 
    # return the field (nphi, ntheta, nr) associated with [varname] 
    # take the jth time slice (default j=0).
    # For "primed" variables, the routine will subtract out the instan-
    # taneous azimuthal average calculated from the equatorial slice
    # _prime in varname refers to az-avg subtracted

    # Get shape of grid
    nphi = eq.nphi
    nr = eq.nr
  
    # Get raw values associated with jth time slice
    vals = eq.vals[:, :, :, j] # now shape (nphi, nr, nq)

    # Get quantity indexes (new or old)---do I actually ever need the 
    # old ones again?
    qind_vr = var_indices['vr']
    qind_vt = var_indices['vt']
    qind_vp = var_indices['vp']
    
    qind_s = var_indices['s']
    qind_p = var_indices['p']
    
    qind_omr = var_indices['omr']
    qind_omt = var_indices['omt']
    qind_omp = var_indices['omp']
    
    qind_br = var_indices['br']
    qind_bt = var_indices['bt']
    qind_bp = var_indices['bp']
    
    if old:
        qind_vr = var_indices_old['vr']
        qind_vt = var_indices_old['vt']
        qind_vp = var_indices_old['vp']
        
        qind_omr = var_indices_old['omr']
        qind_omt = var_indices_old['omt']
        qind_omp = var_indices_old['omp']
        
        qind_s = var_indices_old['s']
        qind_p = var_indices_old['p']       
               
        qind_br = var_indices_old['br']
        qind_bt = var_indices_old['bt']
        qind_bp = var_indices_old['bp']
    
    # Spherical velocities
    if varname in ['vr', 'vr_prime']:
        ind = eq.lut[qind_vr]
        eqslice = vals[:, :, ind]/100. # measure velocities in m/s
    elif varname in ['vt', 'vt_prime']:
        ind = eq.lut[qind_vt]
        eqslice = vals[:, :, ind]/100.
    elif varname in ['vp', 'vp_prime']:
        ind = eq.lut[qind_vp]
        eqslice = vals[:, :, ind]/100.

    # Spherical vorticities
    elif varname in ['omr', 'omr_prime']:
        ind = eq.lut[qind_omr]
        eqslice = vals[:, :, ind]
    elif varname in ['omt', 'omt_prime']:
        ind = eq.lut[qind_omt]
        eqslice = vals[:, :, ind]
    elif varname in ['omp', 'omp_prime']:
        ind = eq.lut[qind_omp]
        eqslice = vals[:, :, ind]
            
    # Enstrophy (total and fluctuating)
    elif varname == 'omsq':
        ind_omr, ind_omt, ind_omp = eq.lut[qind_omr], eq.lut[qind_omt],\
                eq.lut[qind_omp]
        eqslice = vals[:, :, ind_omr]**2 + vals[:, :, ind_omt]**2 +\
                vals[:, :, ind_omp]**2 
    elif varname == 'omsq_prime':
        ind_omr, ind_omt, ind_omp = eq.lut[qind_omr], eq.lut[qind_omt],\
                eq.lut[qind_omp]
        omr_slice = vals[:, :, ind_omr]
        omt_slice = vals[:, :, ind_omt]
        omp_slice = vals[:, :, ind_omp]
        omr_az = np.mean(omr_slice, axis=0).reshape((1, nr))
        omt_az = np.mean(omt_slice, axis=0).reshape((1, nr))
        omp_az = np.mean(omp_slice, axis=0).reshape((1, nr))
        omr_prime = omr_slice - omr_az
        omt_prime = omt_slice - omt_az
        omp_prime = omp_slice - omp_az
        eqslice = omr_prime**2 + omt_prime**2 + omp_prime**2
        
    # Thermodynamic variables: deviations from reference state
    elif varname in ['s', 's_prime']:
        ind = eq.lut[qind_s]
        eqslice = vals[:, :, ind]
    elif varname in ['p', 'p_prime']:
        ind = eq.lut[qind_p]
        eqslice = vals[:, :, ind]

    # Spherical magnetic fields
    elif varname in ['br', 'br_prime']:
        ind = eq.lut[qind_omr]
        eqslice = vals[:, :, ind]
    elif varname in ['bt', 'bt_prime']:
        ind = eq.lut[qind_bt]
        eqslice = vals[:, :, ind]
    elif varname in ['bp', 'bp_prime']:
        ind = eq.lut[qind_bp]
        eqslice = vals[:, :, ind]
    else:
        print("get_eqslice(): unknown variable name %s" %varname)
        print("exiting")
        sys.exit()

    # Subtract out the azimuthal average for the "primed" variables
    if 'prime' in varname and not (varname == 'omsq_prime'):
        eqslice = eqslice - np.mean(eqslice, axis=0).reshape((1, nr))
    return eqslice
