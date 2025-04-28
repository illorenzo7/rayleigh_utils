# Author: Loren Matilsky
# Created: 06/27/2024
#
# Purpose: create a binary file
# (or modify one if it already exists)
# (default name customfile) 
# to contain heating and cooling profiles 
# that are confined to CZ, transitioning to no heating of the RZ 

# Parameters: output_dir (first argument), 

# Command-line options:
#
# --fname: file to (maybe) read reference and save heating
#    (default "customfile")
# --width : width of each the heating and cooling layers:
#   (default 0.10)
# --nr
# Default number of radial (evenly spaced) grid points. 
# Default 10,000 (very fine)

# the following will have defaults set by [fname]_meta.txt, if it exists
# --rmin : bottom of CZ
# --rmax : top of CZ

import numpy as np
import sys, os

sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['raco'] + '/reference_state')

from reference_tools import equation_coefficients
from common import *
from cla_util import *
from HeatingCooling_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas_raw(args)
dirname = clas0['dirname']

# Set default kwargs
kw_default = dotdict(dict({'fname': 'customfile', 'rmin': None, 'rmax': None, 'beta': None, 'nr': None, 'width': 0.1}))
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

# add a heating layer at the bottom and a cooling layer at the top
shape1 = psi_plus(r, kw.rmin, kw.width)
shape2 = psi_minus(r, kw.rmax, kw.width)

# normalize each profile to cancel each other out
A1 = -1. / simps(r**2*shape1, r) # remember rr is in decreasing order
A2 = -1. / simps(r**2*shape2, r)

heat = A1*shape1 - A2*shape2 # this is overall (and smoothed) shape

# now normalize the overall self-cancelling profile

# compute nonradiative flux
Fnr = 1./r**2.*indefinite_integral(heat*r**2., r, kw.rmin)
heat_norm = definite_integral(Fnr*r**2, r, kw.rmin, kw.rmax)
heat_norm /= 1./3.*(kw.rmax**3. - kw.rmin**3.)
heat /= heat_norm

print(buff_line)
print("Computed heating/cooling layers for CZ-only system")
print ("composed of quartics")
print("(rmin, rmax): (%1.5f, %1.5f)" %(kw.rmin, kw.rmax))
print("heating/cooling layer width: %1.5f" %kw.width)
print(buff_line)

# Now write to file using the equation_coefficients framework
print("Setting f_6 in %s" %kw.fname)
# make an executive decision here (since normally f_6 integrates to 1)
# leave f_6 "properly" normalized and later, set c_10 to Ek / Pr or 
# similar
eq.set_function(heat, 6)

print("Writing the heating to %s" %the_file)
print(buff_line)
eq.write(the_file)

# record what we did in the meta file
print("Writing the heating metadata to %s" %meta_file)
print(buff_line)

# we will need the individual lines we want to write
firstline = "Also added custom heating profile using the\n"
lines_new = [firstline]
lines_new.append("generate_HeatingCooling_in_CZ routine.\n")
lines_new.append("heating has the folowing attributes:\n")
lines_new.append("geometry : CZ only\n")
lines_new.append("(rmin, rmax): (%1.5f, %1.5f)\n" %(kw.rmin, kw.rmax))
lines_new.append("width : %1.5f\n" %kw.width)
lastline = buff_line + '\n'
lines_new.append(lastline)
nnew = len(lines_new)

# get the old text
if os.path.isfile(meta_file):
    already_file = True
    f = open(meta_file)
    lines_orig = np.array(f.readlines())
    f.close()
    norig = len(lines_orig)
else:
    already_file = False

# open f again, possibly overwriting heating block (if it exists)
# if not, add the block at the end
f = open(meta_file, 'w')

count = 0
if already_file:
    skip = False
    for line in lines_orig:
        if line == firstline: # we're at the heating block, 
            # start overwriting
            skip = True

        if skip: # write new heating block line
            f.write(lines_new[count])
            count += 1
            if line == lastline:
                skip = False
        else:
            f.write(line)

if not already_file or count == 0: # there was no prior heating block
    for line in lines_new:
        f.write(line)
f.close()
