#/ Author: Loren Matilsky
# Created: 04/26/2024
#
# Purpose: modify a binary file (default name customfile) 
# to contain diffusions (nu, kappa, eta) that go as rho^(power)

# Parameters: output_dir (first argument), 

# Command-line options:
#
# --fname: file to read reference and save heating in (default "customfile")
# --power  : power of rho (defualt -0.5)
# 
# the following will have defaults set by [fname]_meta.txt
# --rmin : bottom of shell
# --rmax : top of shell
# --rt : radius of transition layer

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

# get geometry
jup = False
sun = False
czonly = True

metafile = dirname + '/' + kw_default.fname + '_meta.txt'
f = open(metafile, 'r')
lines = f.readlines()
for line in lines:
    # deal with grid geometry first

    # see if geometry is not CZ only
    alltrue = True
    for keyword in ['rbrz', 'rtrz', 'rbcz', 'rtcz']:
        alltrue *= keyword in line
    if alltrue:
        czonly = False
        sun = True
        st = line.split(':')[1]
        for char in [',', '(', ')']:
            st = st.replace(char, '')
        rmin, rt, rmax = st2 = st.split()

        rmin = float(rmin)
        rt = float(rt)
        rmax = float(rmax)

    if czonly:
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
            rmin = float(rmin)
            rmax = float(rmax)

    if 'Jovian' in line:
        sun = False
        jup = True

# get the top and bottom of CZ
if czonly:
    rbcz, rtcz = rmin, rmax
elif jup:
    rbcz, rtcz = rmin, rt
elif sun:
    rbcz, rtcz = rt, rmin

# add in other default value
kw_default.power = -0.5

# overwrite defaults
kw = update_dict(kw_default, clas)

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# Open and read the hopefully already existing reference file!
eq = equation_coefficients()
the_file = dirname + '/' + kw.fname
eq.read(the_file)
r = eq.radius
nr = eq.nr
rho = eq.functions[0]
dlnrho = eq.functions[7]

# compute diffusion (unnormalized) and its derivative
diffusion = rho**(kw.power)
dlndiffusion = kw.power*dlnrho

# normalize the diffusion
diffusion_norm = definite_integral(diffusion*r**2, r, rbcz, rtcz)
diffusion_norm /= 1./3.*(rtcz**3. - rbcz**3.)
diffusion /= diffusion_norm

print(buff_line)
print("Computed diffusions (nu, kappa, and etea) for RZ-CZ system.")
print("diffusion ~ rho^%0.3f" %kw.power)
print("diffusion normalized by its CZ volume integral")
if jup:
    print ("geometry : Jovian (RZ atop CZ)")
elif sun:
    print ("geometry : solar (CZ atop RZ)")
elif czonly:
    print ("geometry : CZ only")

if czonly:
    print("(rmin, rmax): (%1.5f, %1.5f)" %(rmin, rmax))
else:
    print("(rmin, rt, rmax): (%1.5f, %1.5f, %1.5f)" %(rmin, rt, rmax))
print("power : %1.5f" %kw.power)
print(buff_line)

# Now write to file using the equation_coefficients framework
print("Setting f_3, f_5, and f_7")
for i in [3, 5, 7]:
    eq.set_function(diffusion, i)
print("Setting derivatives f_11, f_12, and f_13")
for i in [11, 12, 13]:
    eq.set_function(dlndiffusion, i)

print("Writing the diffusions to %s" %the_file)
print(buff_line)
eq.write(the_file)

# record what we did in the meta file
f = open(dirname + '/' + metafile, 'a')

f.write("Also added custom diffusion profiles using the\n")
f.write("set_diffusions_rhopower routine.\n")
f.write("diffusions have the folowing attributes:\n")
f.write("diffusion ~ rho^%1.5f\n" %kw.power)
f.write("normalized by CZ volume integral\n")
if jup:
     f.write("geometry : Jovian (RZ atop CZ)\n")
elif sun:
     f.write("geometry : solar (CZ atop RZ)\n")
elif czonly:
     f.write("geometry : CZ only\n")

if czonly:
    f.write("(rmin, rmax): (%1.5f, %1.5f)\n" %(rmin, rmax))
else:
    f.write("(rmin, rt, rmax): (%1.5f, %1.5f, %1.5f)\n" %(rmin, rt, rmax))

f.write("power : %1.5f\n" %kw.power)

f.write(buff_line + '\n')
f.close()
print("Writing the diffusion metadata to %s" %metafile)
print(buff_line)
