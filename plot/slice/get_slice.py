# Author: Loren Matilsky
# Created: 01/03/2019
#
# Extremely long and boring script to find fundamental fluid quantities
# or derivative quantities from a shell slice. Takes a shellslice [a] and
# varname in [vr, vt, vp, vl, vz, ...]
# returns the slice for the variable as an array of shape 
# (nphi, ntheta, nr)
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import *
from varprops import *
from rayleigh_diagnostics import GridInfo # for doing averages

def prime(field): # mean along first axis (phi axis)
    shape = np.shape(field)
    shape_collapsed = np.hstack((np.array([1]), shape[1:]))
    return field - np.mean(field, axis=0).reshape(shape_collapsed)

def prime_sph(field, tw): # doesn't work on equatorial slices
    dummy, nt, nr = np.shape(field)
    field_av = np.mean(field, axis=0) # first take the az-avg
    tw_2d = tw.reshape((nt, 1))
    field_av = np.sum(field_av*tw_2d, axis=0)
    return field - field_av.reshape((1, 1, nr))

def get_slice(a, varname, dirname=None, j=0): 
    # gets a basic field associated with Rayleigh object a
    lut = a.lut 
    rr = a.radius
    # first get the appropriate time slice
    vals = a.vals[..., j]
    if vals.ndim == 4 and not hasattr(a, 'lpower'): 
        # Shell_Slice or Meridional_Slice
        cost = a.costheta
        nt = len(cost)
        cost = cost.reshape((1, nt, 1))
    else:
        # note...don't do cylindrical projections for Shell_Spectra!
        # they don't multiply well
        cost = 0.
        nt = 1
        cost = cost.reshape((nt, nt))

    # first get root variable name and store any modifiers
    varname, deriv, primevar, sphvar = get_varprops(varname)

    # get sine/cotangent from cosine
    sint = np.sin(np.arccos(cost))
    cott = cost/sint
   
    # shape to make geometric fields
    zero = np.zeros(np.array(np.shape(vals[..., 0])))

    # return the basic field based on the variable name
    if varname in var_indices: 
        # this is really only option for Shell_Spectra...
        the_slice = vals[..., lut[var_indices[varname]]]
    elif varname[-1] in ['l', 'z']: # cylindrical variable
        the_slice_r = vals[..., lut[var_indices[varname[:-1] + 'r']]]
        if deriv:
            the_slice_t = vals[..., lut[var_indices[varname[:-1] + 'T']]]
        else:
            the_slice_t = vals[..., lut[var_indices[varname[:-1] + 't']]]
        if varname[-1] == 'l':
            the_slice = sint*the_slice_r + cost*the_slice_t
        elif varname[-1] == 'z':
            the_slice = cost*the_slice_r - sint*the_slice_t
    elif is_an_int(varname):
        the_slice = vals[..., lut[int(varname)]]
    if primevar:
        the_slice = prime(the_slice)
    elif sphvar:
        gi = GridInfo(dirname + '/grid_info')
        tw = gi.tweights
        the_slice = prime_sph(the_slice, tw)
    del vals # free up memory
    return the_slice
