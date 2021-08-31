# Author: Loren Matilsky
# Created: well before 05/06/2019
# Computes the necessary outer entropy gradient to carry out a given
# luminosity (default solar) out of a spherical shell with given 
# reference state and transport coefficients profile
# reference state can either be a polytrope or from a custom file

import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from polytrope import compute_polytrope
from common import *
from cla_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas_raw(args)
dirname = clas0['dirname']

# Set default kwargs
kw_default = dotdict(dict({'lum': lsun, 'rmax': rmax_n3, 'fname': None, 'ktop': 5.0e12, 'rbcz': rbcz, 'rhobcz': rhobcz, 'tempbcz': tempbcz, 'polyn': 1.5, 'polynrho': 3.0, 'cp': c_P}))

# overwrite defaults
kw = update_dict(kw_default, clas)

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# check for common crb kwarg
if 'crb' in clas:
    kw.fname = 'custom_reference_binary'


# get rho*T
if kw.fname is None: # polytrope, default
    print ("got rho and T from POLYTROPIC REFERENCE:")
    print ("poly_nrho : %.2f " %kw.polynrho)
    print ("poly_n    : %.2f " %kw.polyn)
    print ("rho_rbcz  : %.4f " %kw.rhobcz)

    nr = 128 # this is arbitrary -- the top thermodynamic values will be
             # independent of the number of interior grid points
             # We just need arrays (nr > 1) for compute_polytrope to work
             # properly
    di = compute_polytrope(kw.rbcz, kw.rmax, kw.polynrho, nr, kw.polyn, kw.rhobcz)
    rho = di['rho']
    T = di['T']
    rho_rmax = rho[0]
    T_rmax = T[0]
else: # get from custom file
    print('got rho and T from reference file: %s' %kw.fname)
    eq = get_eq(dirname, fname=kw.fname)
    rr = eq.radius
    irmax = np.argmin(np.abs(rr - kw.rmax))
    rho_rmax = eq.rho[irmax]
    T_rmax = eq.T[irmax]
rhot_rmax = rho_rmax*T_rmax 

print ("rho_rmax = %1.8e" %rho_rmax)
print ("T_rmax = %1.8e" %T_rmax)
# Now compute desired entropy gradient
flux_rmax = kw.lum/(4*np.pi*kw.rmax**2)
desired_dsdr = -flux_rmax/rhot_rmax/kw.ktop

print('luminosity : %1.8e erg/s' %kw.lum)
print('kappa_top  : %1.8e cm^2/s' %kw.ktop)
print('rmax       : %1.8e cm' %kw.rmax)
print ('Set outer_dsdr (dtdr_top) to ')
print (make_bold('%1.8e' %desired_dsdr))
