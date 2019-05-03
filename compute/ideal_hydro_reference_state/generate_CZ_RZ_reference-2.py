# Created: 05/01/2019
# Author: Loren Matilsky
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
#

import numpy as np
import sys
from arbitrary_atmosphere import arbitrary_atmosphere
#from struct import pack

import basic_constants as bc

#def pack_array(arr):
#    packed_array = b''
#    n = len(arr)
#    for i in range(n):
#        packed_array += pack('d', arr[i])
#    return packed_array
        
# Set default constants
ri = 3.4139791e10
rm = bc.ri
ro = bc.ro
cp = bc.cp

Tm = bc.T_i
pm = bc.p_i
rhom = bc.rho_i
gam = bc.gamma
k = 1.0
delta = 0.005*rm

# Get directory to save binary files for reference state and heating
dirname = sys.argv[1]

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-ri':
        ri = float(sys.argv[i+1])
    elif arg == '-rm':
        rm = float(sys.argv[i+1])
    elif arg == '-ro':
        ro = float(sys.argv[i+1])
    elif arg == '-nr':
        nr = int(sys.argv[i+1])
    elif arg == '-rhom':
        rhom = float(sys.argv[i+1])        
    elif arg == '-tm':
        tm = float(sys.argv[i+1])   
    elif arg == '-k':
        k = float(sys.argv[i+1])  
    elif arg == '-delta':
        delta = float(sys.argv[i+1])*rm 
    elif arg == '-gam':
        gam = float(sys.argv[i+1]) 
    elif arg == '-cp':
        cp = float(sys.argv[i+1]) 
    elif arg == '-M':
        M = float(sys.argv[i+1])  
        
# First, compute reference state on super-fine grid to interpolate onto later    
nr = 5000
rr = np.linspace(ri, ro, nr)


d2sdr2 = np.zeros_like(rr)
dsdr = k*cp/rm*(1.0 - np.tanh((rr - rm)/delta))
s = k*cp*((rr/rm - 1.0) - (delta/rm)*np.log(np.cosh((rr - rm)/delta)))
g = bc.G*bc.M/rr**2
dgdr = -2.0*g/rr

T, rho, p, dlnT, dlnrho, dlnp, d2lnrho =\
    arbitrary_atmosphere(rr, s, dsdr, d2sdr2, g,\
                         dgdr, rm, Tm, pm, cp, gam)

# Now compute the actual chebyshev colloation points we need
# Look at src/Math_Layer/Chebyshev_Polynomials.F90, lines 165--173 
# (subroutine Initialize_Chebyshev_Grid) for how this works

thefile = dirname + '/custom_reference_binary'
f = open(thefile, "wb")

#may need to specify the data type for a successful read on Rayleigh's end
sigpi = np.array(314,dtype=np.int32)
nr = np.array(nr,dtype=np.int32)
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

#f.write(pack('i', 314))
#f.write(pack('i', nr))
#f.write(pack_array(rr))
#f.write(pack_array(rho))
#f.write(pack_array(dlnrho))
#f.write(pack_array(d2lnrho))
#f.write(pack_array(p))
#f.write(pack_array(T))
#f.write(pack_array(dlnT))
#f.write(pack_array(dsdr))
#f.write(pack_array(s))
#f.write(pack_array(g))
#f.close()