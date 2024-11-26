#/ Author: Loren Matilsky

# Created: 04/26/2024
#
# Purpose: modify a binary file (default name customfile) 
# to contain a ramped up kappa in RZ
# in order to obtain low-sigma regime while being very stiff

# Parameters: output_dir (first argument), 

# Command-line options:
#
# --fname: file to read reference and save heating in (default "customfile")
# --power  : power of rho (default -0.5)
# --delta : thickness of kappa transition (default 0.02)
# --sigma : value of sigma at lower shell boundary (default 0.2)
# --buoy : stiffness value N/2\Omega_0 (default 9.)
# 
# the following will have defaults set by [fname]_meta.txt
# --rmin : bottom of shell
# --rmax : top of shell
# --rt : radius of transition layer
# --jup : if "jup" is specified, RZ lies above CZ

import numpy as np
import sys, os
from generate_CZRZ_reference import psifunc, dpsifunc

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
kw_default = dotdict(dict({'fname': 'customfile', 'delta': 0.02, 'sigma': 0.2, 'buoy': 9.}))
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

# add in other default value kw_default.power = -0.5

# overwrite defaults
kw = update_dict(kw_default, clas)

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# compute geometry of grid
if kw.jup: # RZ above CZ
    rbcz, rtrz = kw.rmin, kw.rmax
    rbrz = rtcz = kw.rt
else: # CZ above RZ
    rbrz, rtcz = kw.rmin, kw.rmax
    rtrz = rbcz = kw.rt

print (rbrz, rtrz, rbcz, rtcz)
if False:
    # Open and read the hopefully already existing reference file!
    eq = equation_coefficients()
    the_file = dirname + '/' + kw.fname
    eq.read(the_file)
    r = eq.radius
    nr = eq.nr
    rho = eq.functions[0]
    dlnrho = eq.functions[7]
    grav = eq.functions[1]/rho
    nsq = grav*eq.functions[13] # make sure I stick with this normalization
    nu = eq.functions[2]
    dlnu = eq.functions[10]
    dnudr = nu*dlnu

    print("nu=", nu)

    if False:
        # get the Prandtl number from main_parameters
        parfile = dirname + '/main_parameters'
        print (buff_line)
        print ("reading Pr from", parfile)
        f = open(parfile, 'r')
        lines = f.readlines()
        f.close()
        di_par = dotdict({})
        for line in lines:
            key, val = line.split('=')
            key = key.strip()
            val = float(val[:-1])
            di_par[key] = val
            print(key, val)
        print("nu=", nu)

        # compute kappa and its derivative to match "sigma" at the lower boundary
        kmax = (nsq*nu)[-1]*kw.buoy*di_par.pr/kw.sigma**2.
        kappa = nu + (kmax - nu[-1])*(1. - psifunc(r, kw.rt - kw.delta, kw.delta))
        numer = dnudr - (kmax - nu[-1])*dpsifunc(r, kw.rt-kw.delta, kw.delta)
        denom = nu + (kmax - nu[-1])*(1. - psifunc(r, kw.rt - kw.delta, kw.delta))
        dlnkappa = numer/denom
        print("nu=", nu)

        # Ok we're done!

        print(buff_line)
        print("Computed diffusion kappa for an RZ-CZ system.")
        print("to ensure arbitrary sigma AND buoyancy via altering Pr(r)")
        if kw.jup:
            print ("geometry : Jovian (RZ atop CZ)")
        else:
            print ("geometry : solar (CZ atop RZ)")
        print("(rmin, rt, rmax): (%1.2f, %1.2f, %1.2f)" %(kw.rmin, kw.rt, kw.rmax))
        print("sigma : %1.3f" %kw.sigma)
        print("buoy : %1.3f" %kw.buoy)
        print("delta_kappa : %1.5f" %kw.delta)
        print(buff_line)

        # Now write to file using the equation_coefficients framework
        print("Setting f_5 and derivative f_12")
        eq.set_function(5, kappa)
        eq.set_function(12, dlnkappa)

        if False:
            print("Writing the diffusions to %s" %the_file)
            print(buff_line)
            eq.write(the_file)

            # record what we did in the meta file
            f = open(dirname + '/' + metafile, 'a')
            print("Also computed a custom diffusion kappa for the RZ-CZ system.")
            print("to ensure arbitrary sigma AND buoyancy via altering Pr(r)")
            if kw.jup:
                print ("geometry : Jovian (RZ atop CZ)")
            else:
                print ("geometry : solar (CZ atop RZ)")
            print("(rmin, rt, rmax): (%1.2f, %1.2f, %1.2f)" %(kw.rmin, kw.rt, kw.rmax))
            print("sigma : %1.3f" %kw.sigma)
            print("buoy : %1.3f" %kw.buoy)
            print("delta_kappa : %1.5f" %kw.delta)
            print(buff_line)


            f.write("Also added custom diffusion profiles using the\n")
            f.write("set_diffusions_rhopower routine.\n")
            f.write("diffusions have the folowing attributes:\n")
            f.write("diffusion ~ rho^%1.5f\n" %kw.power)
            f.write("normalized by CZ volume integral")
            if kw.jup:
                 f.write("geometry : Jovian (RZ atop CZ)\n")
            else:
                 f.write("geometry : solar (CZ atop RZ)\n")
            f.write("(rmin, rt, rmax): (%1.2f, %1.2f, %1.2f)\n" %(kw.rmin, kw.rt, kw.rmax))
            f.write("power : %1.5f\n" %kw.power)

            f.write(buff_line + '\n')
            f.close()
            print("Writing the diffusion metadata to %s" %metafile)
            print(buff_line)
