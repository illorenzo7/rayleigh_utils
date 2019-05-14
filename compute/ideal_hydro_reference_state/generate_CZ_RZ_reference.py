# Author: Loren Matilsky
# Created: 05/01/2019
# Purpose: generate a binary file (for Rayleigh to read) that contains
# a reference state, consisting of a stable region (stiffness k) underlying 
# a neutrally stable CZ

# Parameters: output_dir (first argument), 
# Inner, outer, and transition radii (default ri = 3.4139791e10, 
# rm = 5e10, ro = 6.5860209e10...for CZ corresponding to 3 density
# scale heights )
# nr1 (number of Chebyshevs to use for CZ, -nr1) default 96
# nr2 (number of Chebyshevs to use for RZ, -nr2) default 64
# Density at transition boundary (-rhom) default 0.18053428
# Temperature at transition boundary (-tm) default 2.111256e6
# Stiffness k (-k) default 2
# Transition width delta (-delta) (as a fraction of rm) default 0.005
# gamma (-gam) default 5/3
# pressure specific heat cp (-cp) default 3.5e8
# central mass (-M) default M_sun

import numpy as np
import sys
from arbitrary_atmosphere import arbitrary_atmosphere

import basic_constants as bc

# Set default constants
ri = 4.176e10  # Set RZ width about 0.5x CZ width
rm = bc.ri
ro = bc.ro
cp = bc.cp
rsun = bc.rsun

Tm = bc.T_i
pm = bc.p_i
rhom = bc.rho_i
gam = bc.gamma
k = 2.0
delta = 0.005*rsun

# Get directory to save binary files for reference state and heating
dirname = sys.argv[1]

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-ri':
        ri = float(args[i+1])
    elif arg == '-rm':
        rm = float(args[i+1])
    elif arg == '-ro':
        ro = float(args[i+1])
    elif arg == '-nr':
        nr = int(args[i+1])
    elif arg == '-rhom':
        rhom = float(args[i+1])        
    elif arg == '-tm':
        tm = float(args[i+1])   
    elif arg == '-k':
        k = float(args[i+1])  
    elif arg == '-delta':
        delta = float(args[i+1])*ro
    elif arg == '-gam':
        gam = float(args[i+1]) 
    elif arg == '-cp':
        cp = float(args[i+1]) 
    elif arg == '-M':
        M = float(args[i+1])  
        
# First, compute reference state on super-fine grid to interpolate onto later    
nr = 5000
rr = np.linspace(ri, ro, nr)

# Define a quartic smoothing function:
f = np.zeros(nr)
dfdr = np.zeros(nr)
d2fdr2 = np.zeros(nr)

for i in range(nr):
    rloc = rr[i]
    if rloc <= rm - delta:
        f[i] = 1.0
        dfdr[i] = 0.0
        d2fdr2[i] = 0.0
    elif rloc > rm - delta and rloc < rm:
        f[i] = 1.0 - (1.0 - ((rloc - rm)/delta)**2.0)**3.0
        dfdr[i] = 6.0/delta*(1.0 - ((rloc - rm)/delta)**2.0)**2.0*\
                ((rloc - rm)/delta)
        d2fdr2[i] = 6.0/delta**2.0*(1.0 - ((rloc - rm)/delta)**2.0)*\
                (1.0 - 5.0*((rloc - rm)/delta)**2.0)
    else:
        f[i] = 0.0
        dfdr[i] = 0.0
        d2fdr2[i] = 0.0

line = k*cp*(rr/rm - 1.0)
dlinedr = k*cp/rm*np.ones(nr)
d2linedr2 = np.zeros(nr)

s = line*f
dsdr = dlinedr*f + line*dfdr
d2sdr2 = d2linedr2*f + 2.0*dlinedr*dfdr + line*d2fdr2

g = bc.G*bc.M/rr**2
dgdr = -2.0*g/rr

T, rho, p, dlnT, dlnrho, dlnp, d2lnrho =\
    arbitrary_atmosphere(rr, s, dsdr, d2sdr2, g,\
                         dgdr, rm, Tm, pm, cp, gam)

thefile = dirname + '/custom_reference_binary'
f = open(thefile, "wb")

#may need to specify the data type for a successful read on Rayleigh's end
sigpi = np.array(314, dtype=np.int32)
nr = np.array(nr, dtype=np.int32)
f.write(sigpi.tobytes())
f.write(nr.tobytes())
f.write(rr[::-1].tobytes())
f.write(rho[::-1].tobytes())
f.write(dlnrho[::-1].tobytes())
f.write(d2lnrho[::-1].tobytes())
f.write(p[::-1].tobytes())
f.write(T[::-1].tobytes())
f.write(dlnT[::-1].tobytes())
f.write(dsdr[::-1].tobytes())
f.write(s[::-1].tobytes())
f.write(g[::-1].tobytes())
f.close()
