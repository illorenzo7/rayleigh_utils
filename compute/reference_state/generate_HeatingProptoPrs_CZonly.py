# Author: Loren Matilsky
# Created: 04/26/2024
#
# Purpose: modify a binary file (default name customfile) 
# to contain heating profile that is confined to CZ, transitioning to
# no heating of the RZ -- by default at the stable/unstable layer
# transition
# Must be run AFTER reference state (which includes the 
# density/temperature) is generated

# Parameters: output_dir (first argument), 

# Command-line options:
#
# --fname: file to read reference and save heating in (default "customfile")
# --delta : transition width between CZ and RZ via tanh matching

# the following will have defaults set by [fname]_meta.txt
# --rmin : bottom of shell
# --rmax : top of shell
# --rt : radius of transition layer
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
kw_default = dotdict(dict({'fname': 'customfile', 'rmin': None, 'rmax': None, 'beta': None, 'nr': None}))
# overwrite defaults. Only treat the filename as updated for now, 
# since the other defaults might change
kw = update_dict(kw_default, clas)

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# if there is already a binary file,
# read in its metadata (regarding radial structure) to put into defaults
the_file = dirname + '/' + kw.fname
meta_file = the_file + '_meta.txt'
if os.path.isfile(the_file):
    f = open(meta_file, 'r')
    lines = f.readlines()
    for line in lines:
        # find the line containing rmin, rmax
        alltrue = True
        for keyword in ['rmin', 'rmax']:
            alltrue *= keyword in line
        if alltrue:
            st = line.split(':')[1]
            for char in [',', '(', ')']: # remove parentheses and comma
                st = st.replace(char, '')

            # get rmin and rmax
            rmin, rmax = st.split()

        if 'nr' in line:
            nr = line.split(':')[1]
    f.close()

    # by this point, rmin, rmax, and nr should all be set
    # if not, there was a problem with the meta file
    try:
        kw_default.rmin = float(rmin)
        kw_default.rmax = float(rmax)
        kw_default.nr = float(nr)
    except:
        print ("%s appears to be corrupt. Ignoring its data" %the_file)

# overwrite defaults (again)
kw = update_dict(kw_default, clas)

# if some of the defaults are still None, must overwrite them here
if kw.rmin is None:
    kw.rmin = 3.15
if kw.rmax is None:
    kw.rmax = 4.15

# finally user might have have specified BOTH rmin and rmax via the
# aspect ratio, beta
if not kw.beta is None:
    kw.rmin = kw.beta/(1. - kw.beta)
    kw.rmax = 1./(1. - kw.beta)

# Open and read the possibly already existing reference file
if os.path.isfile(the_file):
    eq = equation_coefficients()
    eq.read(the_file)
    r = eq.radius
    kw.rmin, kw.rmax = np.min(r), np.max(r)
else:
    # compute reference state on super-fine grid to interpolate onto later
    r = np.linspace(kw.rmin, kw.rmax, kw.nr) # keep radius in decreasing order for consistency with Rayleigh convention
    eq = equation_coefficients(r)
nr = len(r)

heat = eq.functions[0]*eq.functions[3] # used to subtract off "top" value in dimensional
    # version. But now realize that top value was at 6.887 x 10^10 cm,
    # way higher than I actually used (so effectively zero). Second
    # paper reported it right. 

# compute nonradiative flux
Fnr = 1./r**2.*indefinite_integral(heat*r**2., r, kw.rmin)
heat_norm = definite_integral(Fnr*r**2, r, kw.rmin, kw.rmax)
heat_norm /= 1./3.*(kw.rmax**3. - kw.rmin**3.)
heat /= heat_norm

print(buff_line)
print("Computed heating for CZ only, profile ~ rho*T")
print("(rmin, rmax): (%1.2f, %1.2f)" %(kw.rmin, kw.rmax))
print(buff_line)

# Now write to file using the equation_coefficients framework
print("Setting f_6")
# make an executive decision here (since normally f_6 integrates to 1)
# leave f_6 "properly" normalized and later, set c_10 to Ek / Pr. 
eq.set_function(heat, 6)

print("Writing the heating to %s" %the_file)
print(buff_line)
eq.write(the_file)

# record what we did in the meta file
f = open(dirname + '/' + meta_file, 'a')

f.write("Also added custom heating profile using the\n")
f.write("generate_CZonly_heating routine.\n")
f.write("Q \propto rho * T.\n")
f.write("heating has the folowing attributes:\n")
f.write("(rmin,  rmax): (%1.2f, %1.2f)\n" %(kw.rmin, kw.rmax))
f.write(buff_line + '\n')
f.close()
print("Writing the heating metadata to %s" %meta_file)
print(buff_line)
