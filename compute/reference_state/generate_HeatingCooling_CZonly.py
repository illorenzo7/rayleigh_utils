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
#   (default 0.05)
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
kw_default.width = 0.05 # this is each heating/cooling layer width
kw_default.nr = 10000

# if there is already a binary file,
# read in its metadata (regarding radial structure) to put into defaults
the_file = dirname + '/' + kw_default.fname
metafile = the_file + '_meta.txt'
if os.path.isfile(the_file):
    f = open(metafile, 'r')
    lines = f.readlines()
    for line in lines:
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
        elif 'solar' in line:
            kw_default.sun = True
            kw_default.rbcz, kw_default.rtcz = rt, rmax
        elif 'Jovian' in line:
            kw_default.jup = True
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
    smooth = 0.5*(1.0 + the_sign*np.tanh((r - kw.rt)/kw.delta)) # "detects" CZ only

# add a heating layer at the bottom and a cooling layer at the top
shape1 = psi_plus(r, rbcz, kw.width)
shape2 = psi_minus(r, rtcz, kw.width)

if kw.jup or kw.sun: # smooth both profiles although it will only matter
    # for the one close to RZ
    shape1 *= smooth
    shape2 *= smooth

# normalize each profile to cancel each other out
A1 = -1. / simps(r**2*shape1, r) # remember rr is in decreasing order
A2 = -1. / simps(r**2*shape2, r)

heat = A1*shape1 - A2*shape2 # this is overall (and smoothed) shape

# now normalize the overall self-cancelling profile

# compute nonradiative flux
Fnr = 1./r**2.*indefinite_integral(heat*r**2., r, rbcz)
heat_norm = definite_integral(Fnr*r**2, r, rbcz, rtcz)
heat_norm /= 1./3.*(rtcz**3. - rbcz**3.)
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
print("Setting f_6 in %s" %fname)
# make an executive decision here (since normally f_6 integrates to 1)
# leave f_6 "properly" normalized and later, set c_10 to Ek / Pr or 
# similar
eq.set_function(heat, 6)

print("Writing the heating to %s" %the_file)
print(buff_line)
eq.write(the_file)

# record what we did in the meta file
print("Writing the heating metadata to %s" %metafile)
print(buff_line)

# get the old text
f = open(dirname + '/' + metafile)
lines_orig = f.readlines()
f.close()

# delete the old custom heating line block
firstline = "Also added custom heating profile using the\n"
skip = False
for line in lines_orig:
    if line == firstline:
        skip = True

    if not skip:
        lines_new.append(line)

    if skip:
        if line == buff_line + '\n':
            skip = False

# re-write the unchanged parts of the meta file
f = open(dirname + '/' + metafile, 'a')

print('lines_new', lines_new)
if False:

    f.write(firstline)
    f.write("generate_HeatingCooling_CZ_only routine.\n")
    f.write("heating has the folowing attributes:\n")
    if kw.jup:
         f.write("geometry : Jovian (RZ atop CZ)\n")
    else:
         f.write("geometry : solar (CZ atop RZ)\n")
    f.write("(rmin, rt, rmax): (%1.2f, %1.2f, %1.2f)\n" %(kw.rmin, kw.rt, kw.rmax))
    f.write("delta_heat : %1.5f\n" %kw.delta)
    f.write("width : %1.5f\n" %kw.width)
    f.write(buff_line + '\n')
    f.close()
