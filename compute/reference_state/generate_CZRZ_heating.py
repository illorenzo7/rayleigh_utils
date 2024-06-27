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
# start with filename, which may change
kw_default = dotdict(dict({'fname': 'customfile'}))
kw_default = update_dict(kw_default, clas)
# read in metadata (regarding radial structure) to put into keywords
metafile = dirname + '/' + kw_default.fname + '_meta.txt'
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

        kw_default.rmin = float(rmin)
        kw_default.rt = float(rt)
        kw_default.rmax = float(rmax)
    elif 'Jovian' in line:
        kw_default.jup = True
f.close()

# add in other default value, the transition width
#kw_default.delta = 0.132
kw_default.delta = 0.005 # make it effectively very sharp, maybe 25 points
# across transition

# overwrite defaults
kw = update_dict(kw_default, clas)

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# compute geometry of grid
if kw.jup: # RZ above CZ
    rbcz, rtrz = kw.rmin, kw.rmax
    rbrz = rtcz = kw.rt
    the_sign = -1.
else: # CZ above RZ
    rbrz, rtcz = kw.rmin, kw.rmax
    rtrz = rbcz = kw.rt
    the_sign = +1.

# Open and read the hopefully already existing reference file!
eq = equation_coefficients()
the_file = dirname + '/' + kw.fname
eq.read(the_file)
r = eq.radius
nr = eq.nr
rho = eq.functions[0]
tmp = eq.functions[3]
smooth = 0.5*(1.0 + the_sign*np.tanh((r - kw.rt)/kw.delta)) # "detects" CZ

heat = rho*tmp*smooth # used to subtract off "top" value in dimensional
    # version. But now realize that top value was at 6.887 x 10^10 cm,
    # way higher than I actually used (so effectively zero). Second
    # paper reported it right. 

# compute nonradiative flux
Fnr = 1./r**2.*indefinite_integral(heat*r**2., r, rbcz)
heat_norm = definite_integral(Fnr*r**2, r, rbcz, rtcz)
heat_norm /= 1./3.*(rtcz**3. - rbcz**3.)
heat /= heat_norm

print(buff_line)
print("Computed heating for RZ-CZ, profile ~ rho*T, joined with a tanh")
if kw.jup:
    print ("geometry : Jovian (RZ atop CZ)")
else:
    print ("geometry : solar (CZ atop RZ)")
print("(rmin, rt, rmax): (%1.2f, %1.2f, %1.2f)" %(kw.rmin, kw.rt, kw.rmax))
print("delta_heat : %1.5f" %kw.delta)
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
f = open(dirname + '/' + metafile, 'a')

f.write("Also added custom heating profile using the\n")
f.write("generate_CZ_RZ_heating-nd routine.\n")
f.write("heating has the folowing attributes:\n")
if kw.jup:
     f.write("geometry : Jovian (RZ atop CZ)\n")
else:
     f.write("geometry : solar (CZ atop RZ)\n")
f.write("(rmin, rt, rmax): (%1.2f, %1.2f, %1.2f)\n" %(kw.rmin, kw.rt, kw.rmax))
f.write("delta_heat : %1.5f\n" %kw.delta)
f.write(buff_line + '\n')
f.close()
print("Writing the heating metadata to %s" %metafile)
print(buff_line)
