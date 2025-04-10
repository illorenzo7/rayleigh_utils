# Author: Loren Matilsky
# Created: 04/24/2024
#
# Purpose: generate a binary file (for Rayleigh to read) that contains
# a reference state, consisting of a stable region
# i.e., a radiative zone (RZ) adjacent to a 
# a neutrally stable CZ

# assumes g(r) = const x 1/r^2

# Parameters: output_dir (first argument), 

# Command-line options:
#
# --beta : ratio of bottom of CZ to top of CZ 
# 
# --gamma : ratio of specific heats
# --nrho : number of density scale heights across CZ
# --fname
# File to save reference state in (default customfile)
#
# --nr
# Default number of radial (evenly spaced) grid points. 
# Default 10,000 (very fine)

import numpy as np
import sys, os
from arbitrary_atmosphere import arbitrary_atmosphere_nd

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
# use double prec for gamma
# 3 sig figs otherwise
kw_default = dotdict(dict({'beta': 0.759, 'gamma': 1.6666666666666667, 'nrho': 3., 'fname': 'customfile', 'nr': 10000}))
# creates profiles from tachocline cases, Matilsky et al. (2022, 2024)

# overwrite defaults
kw = update_dict(kw_default, clas)

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# compute geometry of grid
rmin = kw.beta/(1.-kw.beta)
rmax = 1./(1.-kw.beta)

# note that whatever we get in the end should be rounded to match
# what is manually entered in the main_input file
# for convenience, round radii to two decimal places
rmin = round(rmin, 2)
rmax = round(rmax, 2)

# recompute alpha and beta in light of rounding
kw.beta = rmin/rmax

# compute reference state on super-fine grid to interpolate onto later    
r = np.linspace(rmax, rmin, kw.nr) # keep radius in decreasing order for consistency with Rayleigh convention

# Define an entropy profile that is zero (CZ only)
dsdr = np.zeros(kw.nr)

# compute the atmosphere (this depends on dsdr, not nsq)
rho, tmp, dlnrho, d2lnrho, dlnt, g =\
        arbitrary_atmosphere_nd(r, dsdr, rmin, rmax, kw.gamma, kw.nrho)

# print some info about the atmosphere
print(buff_line)
print("Computed atmosphere for CZ only")
print("nr         : %i" %kw.nr) 
print("beta       : %1.5f" %kw.beta)
print("   (rmin, rmax): (%1.2f, %1.2f)" %(rmin, rmax))
print("gamma      : %1.5f" %kw.gamma)
print("   n=1/(gamma-1)      : %1.5f" %(1./(kw.gamma-1.)))
print("Nrho       : %1.5f" %kw.nrho)
print(buff_line)

# Now write to file using the equation_coefficients framework
eq = equation_coefficients(r)

# Set only the thermodynamic functions/constants in this routine
# In other routines, we can set the heating and transport coefficients
# Only set c_4 = 1/(4*pi) if mag = True

print("Setting f_1, f_2, f_4, f_8, f_9, f_10, and f_14")
eq.set_function(rho, 1)
eq.set_function(rho*g, 2)
eq.set_function(tmp, 4)
eq.set_function(dlnrho, 8)
eq.set_function(d2lnrho, 9)
eq.set_function(dlnt, 10)
eq.set_function(dsdr, 14)

print("Setting all constants to 1.0")
# set all the constants to 1. (no harm I can see for now)
for i in range(eq.nconst):
    eq.set_constant(1.0, i)
the_file = dirname + '/' + kw.fname

print("Writing the atmosphere to %s" %the_file)
eq.write(the_file)

# write metadata to separate file
metafile = the_file + '_meta.txt'
f = open(dirname + '/' + metafile, 'w')

f.write(buff_line + '\n')
f.write("Created custom reference state using the\n")
f.write("generate_CZonly_reference routine.\n")
f.write("The CZ-only system has the folowing attributes:\n")

f.write("nr         : %i\n" %kw.nr) 
f.write("beta       : %1.5f\n" %kw.beta)
f.write("   (rmin, rmax): (%1.2f, %1.2f)" %(rmin, rmax))

# somehow these precisions are confusing me...anyway forget it for now
f.write("gamma      : %1.16f\n" %kw.gamma)
f.write("   n=1/(gamma-1)      : %1.5f\n" %(1./(kw.gamma-1.)))
f.write("Nrho       : %1.5f\n" %kw.nrho)
f.write(buff_line + '\n')
f.close()
print("Writing the metadata   to %s" %metafile)
print(buff_line)
