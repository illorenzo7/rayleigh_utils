# Author: Loren Matilsky
# Created: 05/28/2024
#
# Purpose: estimate (from numerical integrals) the dsdr required
# to balance luminosity from binary file (default name customfile) 
# that contains heating profile that is confined to CZ, transitioning to
# no heating of the RZ 
# Must be run AFTER reference state (which includes the 
# density/temperature/heating/kappa) is generated

# Parameters: output_dir (first argument), 

# Command-line options:
#
# --fname: file to read reference and save heating in (default "customfile")

# the following will have defaults set by [fname]_meta.txt
# --rmin : bottom of shell
# --rmax : top of shell
# --rt : radius of transition layer, if one exists
# --jup : if "jup" is specified, RZ lies above CZ


import numpy as np
import sys, os

sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])

from reference_tools import equation_coefficients
from common import *
from cla_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas_raw(args)
dirname = clas0['dirname']

# Set default kwargs
# start with filename, which may change
kw_default = dotdict(dict({'fname': 'customfile', 'rmin': None, 'rmax': None}))

# overwrite defaults
kw = update_dict(kw_default, clas)

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# Open and read the hopefully already existing reference file!
eq = equation_coefficients()
the_file = dirname + '/' + kw.fname
eq.read(the_file)
r = eq.radius
if kw.rmin is None:
    kw.rmin = np.min(r)
if kw.rmax is None: 
    kw.rmax = np.max(r)
nr = eq.nr
rho = eq.functions[0]
tmp = eq.functions[3]
if eq.fset[4] == 1: # user set kappa via custom ref
    kappa = eq.functions[4]
else: # probably kappa is constant (=1)
    kappa = np.ones_like(rho)
heat = eq.constants[9]*eq.functions[5] # c10*f6
kappa *= eq.constants[5] # kappa is c6*f5

lum = 4.*np.pi*definite_integral(heat*r**2., r, kw.rmin, kw.rmax)
irmax = np.argmin(np.abs(r - kw.rmax))
dtdr_out = -lum/4./np.pi/kw.rmax**2/(rho*tmp*kappa)[irmax]

print ("for lum = int_CZ f_6 dv =  %1.6e" %lum)
print ("where (rmin, rmax) =", kw.rmin, kw.rmax)
print ("set dtdr_top to")
print ("%1.16e" %dtdr_out)
