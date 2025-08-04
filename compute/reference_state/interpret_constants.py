# Author: Loren Matilsky
# Created: 04/27/2024
#
# Purpose: modify binary file with list of human-entered text paramters
# (ra, pr, prm, etc.) and put it into ra_constants list in main_input

# Parameters: output_dir (first argument), 

# Command-line options:
#
# --fname : binary file to modify (default "customfile")
# --tau : timescale in nondimensionalization (default "rot")
#       options rot, visc, kappa
# 

import numpy as np
import sys, os
from arbitrary_atmosphere import compute_Di_v
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
kw_default = dotdict(dict({'fname': 'customfile', 'tau': 'rot'}))

# overwrite defaults
kw = update_dict(kw_default, clas)

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# get the name of the metafile
metafile = dirname + '/' + kw.fname + '_meta.txt'

# get paremeters necessary for dissipation number
f = open(metafile, 'r')
lines = f.readlines()

# deal with geometry
jup = False
sun = False
czonly = True
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
    elif 'gamma' in line and not 'n' in line: # gamma appears twice...
        gamma = float(line.split(':')[1][:-1])
    elif 'Nrho' in line:
        nrho = float(line.split(':')[1][:-1])
f.close()

# recompute beta in light of rounding
if jup: # RZ above CZ
    beta = rmin/rt
elif sun: # CZ above RZ
    beta = rt/rmax
elif czonly:
    beta = rmin/rmax

# now get dissipation number
di = compute_Di_v(gamma, beta, nrho)

# read in parameters from main_parameters
parfile = dirname + '/main_parameters'
print (buff_line)
print ("reading", parfile)
f = open(parfile, 'r')
lines = f.readlines()
f.close()
di_par = dotdict({})
for line in lines:
    key, val = line.split('=')
    key = key.strip()
    val = float(val[:-1])
    di_par[key] = val
    print (key, "=", val)

print ("(gamma, beta, Nrho) = (%1.16f, %1.5f, %1.5f)" %(gamma, beta, nrho))
print ("   ----> Di = %1.5f" %di)

# "unpack" pr and ra (or supercrit)
pr = di_par.pr

if 'supercrit' in di_par.keys():
    supercrit = di_par.supercrit
    rotation = True
    raf = supercrit**3/pr**2/di_par.roc**4
    ek = di_par.roc * np.sqrt(pr/raf)
else:
    raf = di_par.raf

    # if model is rotating, user will have set Ro_c
    # must compute ekman number in that case
    if 'roc' in di_par.keys():
        rotation = True
        ek = di_par.roc * np.sqrt(pr/raf)
        rafmod = di_par.roc**2.
        supercrit = raf*ek**(4/3)
    else:
        rotation = False

# if model has stable layer 
# user will have set sigma or bu
# must compute buoyancy parameter in that case
if 'sigma' in di_par.keys(): # then the model is rotating
    bumod = di_par.sigma**2./pr
    bu = bumod*pr/ek**2.
elif 'bu' in di_par.keys():
    bu = di_par.bu
else: # no stable layer
    bu = 0.

if 'prm' in di_par.keys():
    magnetism = True # this is just to set the constants
    # the same reference file can be used for a hydro run
    prm = di_par.prm
else:
    magnetism = False

# compute tau over taunu [will have to update if I change which parameter
# set I choose by default...
# right now its Ro_c, sigma, Ra, Pr, Pr_m for rotating cases,
# etc.
if kw.tau == 'rot':
    timescale_msg = "timescale chosen: %s, tau = 1/(2 omega_0)" %kw.tau
    tau_over_taunu = ek
    print ("convective Rossby number = %1.5e" %di_par.roc)
    print ("Ekman number = %1.5e" %ek)
    print ("Rayleigh number = %1.5e" %raf)
    print ("supercriticality = %1.5e" %supercrit)
elif kw.tau == 'visc':
    timescale_msg ="timescale chosen: %s; tau = H^2/nu" %kw.tau
    tau_over_taunu = 1.
elif kw.tau == 'kappa':
    timescale_msg = "timescale chosen: %s; tau = H^2/kappa" %kw.tau
    tau_over_taunu = pr
print(timescale_msg)

# read the existing equation coefficients file
eq = equation_coefficients()
the_file = dirname + '/' + kw.fname
eq.read(the_file)

# now set all the constants

# c1
if rotation:
    eq.set_constant(tau_over_taunu/ek, 1)
else:
    eq.set_constant(0., 1)

# c2
eq.set_constant(raf/pr*tau_over_taunu**2., 2)

# c3
eq.set_constant(1., 3)

# c4
if magnetism:
    eq.set_constant(1., 4)
else:
    eq.set_constant(0., 4)

# c5
eq.set_constant(tau_over_taunu, 5)

# c6
eq.set_constant(tau_over_taunu/pr, 6)

# c7
if magnetism:
    eq.set_constant(tau_over_taunu/prm, 7)
else:
    eq.set_constant(0., 7)

# c8
eq.set_constant(pr*di/raf/tau_over_taunu**2., 8)

# c9
if magnetism:
    eq.set_constant(pr*di/raf/tau_over_taunu**2., 9)
else:
    eq.set_constant(0., 9)

# c10
eq.set_constant(tau_over_taunu/pr, 10)

# c11
eq.set_constant(bu/raf, 11)

# write the file
print("Writing the constants to %s" %the_file)
eq.write(the_file)

# record what we did in the meta file
print("and saying I did so in the metadata file %s" %metafile)
print(buff_line)

f = open(dirname + '/' + metafile, 'a')

f.write("Also set constants c1 thru c11 using\n")
f.write("%s ---> %s\n" %(parfile, the_file))
f.write("with the interpret_constants routine.\n")
f.write(timescale_msg + '\n')
f.write("Di: %1.5f\n" %di)
if kw.tau == 'rot':
    f.write("convective Rossby number: %1.5e\n" %di_par.roc)
    f.write("Ekman number: %1.5e\n" %ek)
    f.write("Rayleigh number = %1.5e\n" %raf)
    f.write("supercritcality: %1.5e\n" %supercrit)

f.write(buff_line + '\n')
f.close()
