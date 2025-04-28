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
# --delta : possible transition width between CZ and RZ via tanh matching
#   (default quite sharp: 0.005)
# --width : width of each the heating and cooling layers:
#   (default 0.10)
# --nr
# Default number of radial (evenly spaced) grid points. 
# Default 10,000 (very fine)

# the following will have defaults set by [fname]_meta.txt, if it exists
# --rbcz : bottom of CZ
# --rtcz : top of CZ
# --jup : (default False or what's in customfile) 
#        if "jup" is True, there is RZ above CZ, use smoothing
# --sun : (default False or what's in customfile) 
#        if "sun" is True, there is RZ below CZ, use smoothing

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
# start with filename, which may change
kw_default = dotdict(dict({'fname': 'customfile', 'jup': False, 'sun': False, 'rbcz': None, 'rtcz': None}))

# add in other default value, the transition width from CZ to RZ
kw_default.delta = 0.005 # make it effectively very sharp, maybe 25 points
kw_default.width = 0.1 # this is each heating/cooling layer width
kw_default.nr = 10000

# if there is already a binary file,
# read in its metadata (regarding radial structure) to put into defaults
the_file = dirname + '/' + kw_default.fname
meta_file = the_file + '_meta.txt'
if os.path.isfile(the_file):
    f = open(meta_file, 'r')
    lines = f.readlines()
    for line in lines:
        # see if geometry might be solar-like
        if 'solar' in line:
            kw_default.sun = True
        elif 'Jovian' in line:
            kw_default.jup = True

        # find the line containing rbrz, etc.
        alltrue = True
        for keyword in ['rbrz', 'rtrz', 'rbcz', 'rtcz']:
            alltrue *= keyword in line
        if alltrue:
            st = line.split(':')[1]
            for char in [',', '(', ')']:
                st = st.replace(char, '')
            rmin, rt, rmax = st2 = st.split()

            rmin = float(rmin)
            rt = float(rt)
            rmax = float(rmax)
            if kw_default.sun:
                kw_default.rbcz, kw_default.rtcz = rt, rmax
            if kw_default.jup:
                kw_default.rbcz, kw_default.rtcz = rmin, rt
    f.close()

# overwrite defaults
kw = update_dict(kw_default, clas)

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# Open and read the possibly already existing reference file
if os.path.isfile(the_file):
    eq = equation_coefficients()
    eq.read(the_file)
    r = eq.radius
else:
    # compute reference state on super-fine grid to interpolate onto later
    r = np.linspace(kw.rtcz, kw.rtrz, kw.nr) # keep radius in decreasing order for consistency with Rayleigh convention
    # also, this logic assumes that no prior file means "CZ only"
    eq = equation_coefficients(r)
nr = len(r)

# if there is an RZ, compute a smoothing profile
if kw.jup: # there is RZ above CZ
    rt = kw.rtcz
    the_sign = -1.
elif kw.sun: # there is RZ below CZ
    rt = kw.rbcz
    the_sign = +1.
if kw.jup or kw.sun:
    smooth = 0.5*(1.0 + the_sign*np.tanh((r - rt)/kw.delta)) # "detects" CZ only

# add a heating layer at the bottom and a cooling layer at the top
shape1 = psi_plus(r, kw.rbcz, kw.width)
shape2 = psi_minus(r, kw.rtcz, kw.width)

if kw.jup or kw.sun: # smooth both profiles although it will only matter
    # for the one close to RZ
    shape1 *= smooth
    shape2 *= smooth

# normalize each profile to cancel each other out
A1 = -1. / simpson(r**2*shape1, x=r) # remember rr is in decreasing order
A2 = -1. / simpson(r**2*shape2, x=r)

heat = A1*shape1 - A2*shape2 # this is overall (and smoothed) shape

# now normalize the overall self-cancelling profile

# compute nonradiative flux
Fnr = 1./r**2.*indefinite_integral(heat*r**2., r, kw.rbcz)
heat_norm = definite_integral(Fnr*r**2, r, kw.rbcz, kw.rtcz)
heat_norm /= 1./3.*(kw.rtcz**3. - kw.rbcz**3.)
heat /= heat_norm

print(buff_line)
print("Computed heating/cooling layers in CZ, made with quartics")
if kw.jup:
    print ("geometry : Jovian (RZ atop CZ)")
elif kw.sun:
    print ("geometry : solar (CZ atop RZ)")
else:
    print ("geometry : CZ only")
if kw.jup or kw.sun:
    print("(rmin, rt, rmax): (%1.2f, %1.2f, %1.2f)" %(rmin, rt, rmax))
else:
    print("(rmin, rmax): (%1.2f, %1.2f)" %(kw.rbcz, kw.rtcz))
print("delta_heat : %1.5f" %kw.delta)
print("heating/cooling layer width_heat : %1.5f" %kw.width)
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
if kw.jup:
     lines_new.append("geometry : Jovian (RZ atop CZ)\n")
elif kw.sun:
     lines_new.append("geometry : solar (CZ atop RZ)\n")
else:
     lines_new.append("geometry : CZ only\n")
if kw.jup or kw.sun:
    lines_new.append("(rmin, rt, rmax): (%1.2f, %1.2f, %1.2f)\n" %(rmin, rt, rmax))
else:
    lines_new.append("(rmin, rmax): (%1.2f, %1.2f)\n" %(rmin, rmax))
lines_new.append("delta_heat : %1.5f\n" %kw.delta)
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
