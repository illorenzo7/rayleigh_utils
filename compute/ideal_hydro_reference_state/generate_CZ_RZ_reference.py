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
from scipy.interpolate import BSpline, splrep
from struct import pack

import basic_constants as bc

def spline(xk, yk, xnew):
    # Python in its infinite wisdom deprecated the nice "spline" routine
    # in scipy.interpolate--this new routine mimics the old "spline" syntax
    # but uses BSpline classes
    t, c, k = splrep(xk, yk)
    spline_loc = BSpline(t, c, k)
    return spline_loc(xnew)

def pack_array(arr):
    packed_array = b''
    n = len(arr)
    for i in range(n):
        packed_array += pack('d', arr[i])
    return packed_array
        
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
nr1 = 96    # domain 1 == CZ
nr2 = 64    # domain 1 == RZ

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
r_fine = np.linspace(ri, ro, 1000)


d2sdr2_fine = np.zeros_like(r_fine)
dsdr_fine = k*cp/rm*(1.0 - np.tanh((r_fine - rm)/delta))
s_fine = k*cp*((r_fine/rm - 1.0) - (delta/rm)*np.log(np.cosh((r_fine - rm)/delta)))
g_fine = bc.G*bc.M/r_fine**2
dgdr_fine = -2.0*g_fine/r_fine

T_fine, rho_fine, p_fine, dlnT_fine, dlnrho_fine, dlnp_fine, d2lnrho_fine =\
    arbitrary_atmosphere(r_fine, s_fine, dsdr_fine, d2sdr2_fine, g_fine,\
                         dgdr_fine, rm, Tm, pm, cp, gam)

# Now compute the actual chebyshev colloation points we need
# Look at src/Math_Layer/Chebyshev_Polynomials.F90, lines 165--173 
# (subroutine Initialize_Chebyshev_Grid) for how this works
    
dth = np.pi/nr1
theta = np.linspace(dth/2, np.pi - dth/2, nr1)
x = np.cos(theta)
scaling = (ro - rm)/(x[0] - x[-1])
r1 = rm + (x - x[-1])*scaling

dth = np.pi/nr2
theta = np.linspace(dth/2, np.pi - dth/2, nr2)
x = np.cos(theta)
scaling = (rm - ri)/(x[0] - x[-1])
r2 = ri + (x - x[-1])*scaling

rr = np.hstack((r1, r2))

T = spline(r_fine, T_fine, rr)
rho = spline(r_fine, rho_fine, rr)
p = spline(r_fine, p_fine, rr)
dlnT = spline(r_fine, dlnT_fine, rr)
dlnrho = spline(r_fine, dlnrho_fine, rr)
d2lnrho = spline(r_fine, d2lnrho_fine, rr)
s = spline(r_fine, s_fine, rr)
dsdr = spline(r_fine, dsdr_fine, rr)
g = spline(r_fine, g_fine, rr)

thefile = dirname + '/custom_reference_binary'
f = open(thefile, "wb")
sigpi = np.array(314, dtype=np.int32)
nr = np.array(nr1 + nr2, dtype=np.int32)
f.write(pack('i', 314))
f.write(pack('i', nr1 + nr2))
f.write(pack_array(rr))
f.write(pack_array(rho))
f.write(pack_array(dlnrho))
f.write(pack_array(d2lnrho))
f.write(pack_array(p))
f.write(pack_array(T))
f.write(pack_array(dlnT))
f.write(pack_array(dsdr))
f.write(pack_array(s))
f.write(pack_array(g))
f.close()