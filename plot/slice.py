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
from rayleigh_diagnostics import GridInfo # for doing averages

var_indices = {\
    'vr'    :       1, 
    'vt'    :       2,
    'vp'    :       3,              

    'dvrdr' :       10,
    'dvtdr' :       11,
    'dvpdr' :       12,

    'dvrdt' :       19,
    'dvtdt' :       20,
    'dvpdt' :       21,

    'dvrdp' :       28,
    'dvtdp' :       29,
    'dvpdp' :       30,

    'dvrdT' :       37,
    'dvtdT' :       38,
    'dvpdT' :       39,

    'dvrdP' :       46,
    'dvtdP' :       47,
    'dvpdP' :       48,

    'omr'   :       301,
    'omt'   :       302,
    'omp'   :       303,

    's'     :       501,
    'p'     :       502,

    'br'    :       801,
    'bt'    :       802,
    'bp'    :       803,

    'br'    :       1, 
    'bt'    :       2,
    'bp'    :       3,              

    'dbrdr' :       810,
    'dbtdr' :       811,
    'dbpdr' :       812,

    'dbrdt' :       819,
    'dbtdt' :       820,
    'dbpdt' :       821,

    'dbrdp' :       828,
    'dbtdp' :       829,
    'dbpdp' :       830,

    'dbrdT' :       837,
    'dbtdT' :       838,
    'dbpdT' :       839,

    'dbrdP' :       846,
    'dbtdP' :       847,
    'dbpdP' :       848,

    'jr'    :       1001,
    'jt'    :       1002,
    'jp'    :       1003}

rootlabels = {'v': r'$v$', 'b': r'$B$', 'om': r'$\omega$', 'j': r'$\mathcal{J}$'}
dirlabels = {'r': r'$_r$', 't': r'$_\theta$', 'p': r'$_\phi$', 'l': r'$_\lambda$', 'z': r'$_z$', 'T': r'$_\Theta$', 'P': r'$_\Phi$'}

def get_varprops(varname):
    # get boolean values deriv, prime (az average subtracted), and sph
    # (sph average subtracted)
    deriv = False
    prime = False
    sph = False
    if 'dr' in varname or 'dt' in varname or 'dp' in varname or 'dT' in varname or 'dP' in varname:
        deriv = True
    if varname[-6:] == '_prime': # prime appears at end to modify varname
        prime = True
        varname = varname[:-6]
    elif varname[-4:] == '_sph'
        sph = True
        varname = varname[:-4]
    return varname, deriv, prime, sph

def get_fieldlabel(varname):
    varname, deriv, prime, sph = compute_varprops(varname)
    if deriv:
        derivdir = varname[-1]
        varname = varname[1:]
        varname = varname[:-2]
    # now get root label
    # start with thermal vars
    if varname == 'p':
        label = r'$P$'
    elif varname == 's':
        label = r'$S$'
    elif varname == 't':
        label = r'$T$'
    elif varname == 'rho':
        label = r'$\rho$'
    elif varname == 'cost':
        label = r'$\cos\theta$'
    elif varname == 'sint':
        label = r'$\sin\theta$'
    elif varname == 'cott':
        label = r'$\cot\theta$'
    elif varname == 'rr':
        label = r'$r$'
    # now field variables, with a vector component
    if 'v' in varname or 'b' in varname or 'om' in varname or 'j' in varname:
        rootname = varname[:-1]
        direction = varname[-1]
        label = rootlabels[rootname] + dirlabels[direction]
    if prime:
        label += r'$^\prime$'
    elif sph:
        label += r'$^{\prime\prime}$'
    if deriv:
        label = r'$\partial$' + label + r'$/$' + r'$\partial$' + dirlabels[derivdir]
    return label

def prime(field): # mean along first axis (phi axis)
    shape = np.shape(field)
    shape_collapsed = np.hstack((np.array([1]),shape[1:]))
    return field - np.mean(field, axis=0).reshape(shape_collapsed)

def prime_sph(field, tw): # doesn't work on equatorial slices
    dummy, nt, nr = np.shape(field)
    field_av = np.mean(field, axis=0) # first take the az-avg
    tw_2d = tw.reshape((nt, 1))
    field_av = np.sum(field_av*tw_2d, axis=0)
    return field - field_av.reshape((1, 1, nr))

var_indices = {\
    'vr'    :       1, 
    'vt'    :       2,
    'vp'    :       3,              

    'dvrdr' :       10,
    'dvtdr' :       11,
    'dvpdr' :       12,

    'dvrdt' :       19,
    'dvtdt' :       20,
    'dvpdt' :       21,

    'dvrdp' :       28,
    'dvtdp' :       29,
    'dvpdp' :       30,

    'dvrdT' :       37,
    'dvtdT' :       38,
    'dvpdT' :       39,

    'dvrdP' :       46,
    'dvtdP' :       47,
    'dvpdP' :       48,

    'omr'   :       301,
    'omt'   :       302,
    'omp'   :       303,

    's'     :       501,
    'p'     :       502,

    'br'    :       801,
    'bt'    :       802,
    'bp'    :       803,

    'br'    :       1, 
    'bt'    :       2,
    'bp'    :       3,              

    'dbrdr' :       810,
    'dbtdr' :       811,
    'dbpdr' :       812,

    'dbrdt' :       819,
    'dbtdt' :       820,
    'dbpdt' :       821,

    'dbrdp' :       828,
    'dbtdp' :       829,
    'dbpdp' :       830,

    'dbrdT' :       837,
    'dbtdT' :       838,
    'dbpdT' :       839,

    'dbrdP' :       846,
    'dbtdP' :       847,
    'dbpdP' :       848,

    'jr'    :       1001,
    'jt'    :       1002,
    'jp'    :       1003}

rootlabels = {'v': r'$v$', 'b': r'$B$', 'om': r'$\omega$', 'j': r'$\mathcal{J}$'}
dirlabels = {'r': r'$_r$', 't': r'$_\theta$', 'p': r'$_\phi$', 'l': r'$_\lambda$', 'z': r'$_z$', 'T': r'$_\Theta$', 'P': r'$_\Phi$'}
operator_symbols = ['+', '-', '*', '/', '(', ')']

def compute_texlabel(varname):
    deriv = False
    prime = False
    sph = False

    if 'dr' in varname or 'dt' in varname or 'dp' in varname or 'dT' in varname or 'dP' in varname:
        deriv = True
        derivdir = varname[-1]
        varname = varname[1:]
        varname = varname[:-2]

    if 'prime' in varname: # prime appears at end to modify varname
        prime = True
        varname = varname[:-5]
    elif 'sph' in varname:
        sph = True
        varname = varname[:-3]

    # now get root label
    # start with thermal vars
    if varname == 'p':
        label = r'$P$'
    elif varname == 's':
        label = r'$S$'
    elif varname == 't':
        label = r'$T$'
    elif varname == 'rho':
        label = r'$\rho$'
    
    # now field variables, with a vector component
    if 'v' in varname or 'b' in varname or 'om' in varname or 'j' in varname:
        rootname = varname[:-1]
        direction = varname[-1]
        label = rootlabels[rootname] + dirlabels[direction]

    if prime:
        label += r'$^\prime$'
    elif sph:
        label += r'$^{\prime\prime}$'
    if deriv:
        label = r'$\partial$' + label + r'$/$' + r'$\partial$' + dirlabels[derivdir]
    return label

def smooth(field, dlon):
    nphi, nt, nr = np.shape(field)
    nphi_av = int(nphi*dlon/360.) + 1 # averaging over this number will
        # yield an interval as close as possible to dlon
    ov2 = nphi_av//2
    field_smooth = np.zeros_like(field)
    # calculate the smoothed field phi-index by phi-index
    iphimin = -ov2
    if nphi_av % 2 == 0: # nphi_av is even
        iphimax = ov2
    else: # nphi_av is odd
        iphimax = ov2 + 1
    for i in np.arange(iphimin, iphimax):
        field_smooth[0] += field[i]
    for i in range(1, nphi):
        to_sub = field[iphimin + i - 1]
        to_add = field[(iphimax + i - 1)%nphi]
        field_smooth[i] += (field_smooth[i-1] + to_add - to_sub)
    return field_smooth/nphi_av

def get_field(vals, lut, rr, cost, varname, dirname=None): 
    # gets the basic field associated with vals array
    # more complicated slices can be achieved by doing arithmetic
    # operations on the individual field
    # get geometric terms
    # if this is one of the standard field variables we're done!

    # first get root variable name and store any modifiers
    varname, deriv, prime, sph = compute_varprops(varname)

    # get sine/cotangent from cosine
    sint = np.sin(np.acos(cost))
    cott = cost/sint
   
    # shape to make geometric fields
    zero = np.zeros(np.array(np.shape(vals[..., 0])))

    # return the basic field based on the variable name
    if varname in var_indices.keys():
        the_field = vals[..., lut[var_indices[varname]]]
    elif varname == 'cost':
        the_field = zero + cost
    elif varname == 'sint':
        the_field = zero + sint
    elif varname == 'cott':
        the_field = zero + cott
    elif varname == 'rr':
        the_field = zero + rr
    elif varname in ['rho', 't']: # derived thermo variable
        # need background thermo reference
        eq = get_eq(dirname)
        rr_full = eq.radius # differs from rr for Shell_Slices
        shape_collapsed_r = np.shape(rr)
        rr_1d = rr[..., :]
        nr = len(rr_1d)
        rinds = np.zeros(nr, dtype='int')
        for ir in range(nr):
            rinds[ir] = np.argmin(np.abs(rr_full - rr[ir]))
        ref_rho = (eq.density)[rinds].reshape(shape_collapsed_r)
        ref_T = (eq.temperature)[rinds].reshape(shape_collapsed_r)
        ref_P = (eq.pressure)[rinds].reshape(shape_collapsed_r)
        s_field = vals[..., lut[var_indices['s']]]
        p_field = vals[..., lut[var_indices['p']]]
        if varname == 'rho':
            the_field = ref_rho*(p_field/ref_P/thermo_gamma - s_field/c_P) 
        elif varname == 't':
            the_field = ref_T*(p_field/ref_P*(1. - 1./thermo_gamma) +\
                    s_field/c_P)
    elif varname[-1] in ['l', 'z']: # cylindrical variable
        the_field_r = vals[..., lut[var_indices[varname[:-1] + 'r']]]
        if deriv:
            the_field_t = vals[..., lut[var_indices[varname[:-1] + 'T']]]
        else:
            the_field_t = vals[..., lut[var_indices[varname[:-1] + 't']]]
        if varname[-1] == 'l':
            the_field = sint*the_field_r + cost*the_field_t
        elif varname[-1] == 'z':
            the_field = cost*the_field_r - sint*the_field_t
    elif is_an_int(varname):
        the_field = vals[..., lut[int(varname)]]
    if prime:
        the_field = prime(the_field)
    elif sph:
        gi = GridInfo(dirname + '/grid_info')
        tw = gi.tweights
        the_field = prime_sph(the_field, tw)
    del vals # free up memory
    return the_field

def resolve_expression(vals, lut, rr, cost, expr, dirname=None, starting_expression=None):
    shape = np.array(np.shape(vals[..., 0]))
    field1 = np.zeros(shape)
    for st1 in expr.split('+'):
        field2 = np.zeros(shape)
        for st2 in st1.split('-'):
            field3 = np.ones(shape)
            for st3 in st2.split('*'):
                count = 1
                for varname in st3.split('/'):
                    if varname == 'expr':
                        field4 = starting_expression
                    else:
                        if count == 1:
                            field4 = get_field(vals, lut, rr, cost, varname, dirname=dirname)
                        else:
                            field4 /= get_field(vals, lut, rr, cost, varname, dirname=dirname)
                    count += 1
                field3 *= field4
            field2 -= field3
        field1 += field2
    del vals # free up memory
    return field1 

    # possibly smooth the variable
    if smooth_desired:
        sslice = smooth(sslice, nphi_av)

def get_slice(a, varname, dirname=None, j=0):
    # Given a slice object (vals array shape nphi, ntheta, nr, nq, niter=0)
    # return the slice (nphi, ntheta, nr) associated with [varname]
    # can do arithmetic operations on the basic fields

    # first get root variable name and store any modifiers
    varname, deriv, prime, sph = compute_varprops(varname)

    # Get raw values associated with jth time slice
    vals = a.vals[..., j]
    lut = a.lut
    
    # grid_info
    shape = np.array(np.shape(vals[..., 0]))
    if len(shape) == 3:
        nphi, nt, nr = shape
        shape_collapsed = np.array([1, nt, 1])
        cost = a.sintheta.reshape(shape_collapsed)
        sint = a.sintheta.reshape(shape_collapsed)
        cott = cost/sint
        shape_collapsed_r = np.array([1, 1, nr])
    elif len(shape) == 2: # Equatorial slice
        nphi, nr = shape
        cost = 0.
        sint = 1.
        cott = 0.
        shape_collapsed_r = np.array([1, 1, nr])
    rr = a.radius.reshape(shape_collapsed_r)

    # see if we will smooth over longitude at the end
    smooth_desired = False
    if 'smooth' in varname:
        smooth_desired = True
        varname = varname.replace('smooth', '')
        nphi_av = int(varname[:3])
        varname = varname[3:]

    # make sure expression is bracketed by parentheses
    varname = '(' + varname + ')'

    # resolve expressions in parentheses, starting from innermost set
    lenvar = len(varname)

    # Get locations of parentheses
    ileft = []
    iright = []
    nparenth = 0 # number of parentheses
    # get left ones first
    for i in range(lenvar):
        if varname[i] == '(':
            nparenth += 1
            ileft.append(i)
    # now right ones; remember leftmost left corresponds to rightmost right
    nparenth_right = 0
    for i in range(lenvar - 1, -1, -1):
        if varname[i] == ')':
            nparenth_right += 1
            iright.append(i)
    if nparenth_right != nparenth: # someone done fucked up!
        print ("Invalid expression %i left parenth., but %i right parenth."\
                %(nparenth, nparenth_right))
        sys.exit()
   
    # do arithmetic; resolve the innermost parentheses first
    for i in range(nparenth - 1, -1, -1):
        ileft_loc = ileft[i]
        iright_loc = iright[i]
        expr = varname[ileft_loc + 1: iright_loc]
        if i == nparenth - 1: # innermost expression
            the_slice = resolve_expression(vals, lut, rr, cost, expr, dirname=None, starting_expression=None):
        else:
            expr = varname.replace('(' + old_expr + ')', 'expr')
            expr = varname_sh[ileft_loc + 1:iright_loc]
            the_slice = resolve_expression(vals, lut, rr, cost, expr.replace(old_expr, 'expr'), dirname=None, starting_expression=the_slice):
        old_expr = expr
    if prime:
        the_slice = prime(the_slice)
    elif sph:
        gi = GridInfo(dirname + '/grid_info')
        tw = gi.tweights
        the_slice = prime_sph(the_slice, tw)
    del vals # free up memory
    return the_slice
