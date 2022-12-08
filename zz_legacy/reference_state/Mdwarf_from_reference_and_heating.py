# Author: Loren Matilsky
# Created: 11/26/2019
#
# Purpose: generate a binary file (for Rayleigh to read) that contains
# a reference state and heating function, from one of Connor's old
# custom reference states

# Parameters: 
# output_dir (first argument), 
# reference_file (usually mesa_reference)
# heating_file (usually mesa_heating)

import numpy as np
import sys, os
import struct
from scipy.integrate import cumtrapz

sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['raco'] + '/tachocline')

from reference_tools import equation_coefficients
from rayleigh_diagnostics import ReferenceState

cp = 3.5e8

# Get directory to save binary files for reference state and heating
dirname = sys.argv[1]
reference_file = sys.argv[2]
heating_file = sys.argv[3]

fname = 'custom_reference_binary'

ref = ReferenceState(reference_file, '')

nr = ref.nr
rr = ref.radius

print ("Read reference stuff from ", reference_file)

heating_data = open(heating_file, 'rb').read()
fmt_st = '%id' %nr
heating = struct.unpack(fmt_st, heating_data[8+nr*8:8+nr*8*2])
rhot = ref.density*ref.temperature
rmin, rmax = 9.07e9, 2.507e10
ir1, ir2 = np.argmin(np.abs(rr - rmax)), np.argmin(np.abs(rr - rmin))
heating_int = -4.*np.pi*cumtrapz((rr**2*rhot*heating)[ir1:ir2+1],\
        rr[ir1:ir2+1])[-1]
print('heating int is ', heating_int)
f6 = rhot*heating/heating_int

print("read heating data from ", heating_file)

ldwarf = 9.478e31

# Now write to file using the equation_coefficients framework
eq = equation_coefficients(rr)

# Set only the thermodynamic functions/constants in this routine
# In other routines, we can set the heating and transport coefficients
# Only set c_4 = 1/(4*pi) if mag = True

print("Setting f_1, f_2, f_4, f_8, f_9, f_10, and f_14")
eq.set_function(ref.density, 1)
buoy = ref.density*ref.gravity/cp
eq.set_function(buoy, 2)
eq.set_function(ref.temperature, 4)
eq.set_function(ref.dlnrho, 8)
eq.set_function(ref.d2lnrho, 9)
eq.set_function(ref.dlnt, 10)
eq.set_function(ref.dsdr, 14)

print("Setting c_2 and c_3")
print("magnetism = False, so setting c_4 = 0.0")
eq.set_constant(0.0, 4) # multiplies Lorentz force

eq.set_constant(1.0, 2) # multiplies buoyancy
eq.set_constant(1.0, 3) # multiplies pressure grad.

print ("Setting heating: c_10 and f_6")
eq.set_function(f6, 6)
eq.set_constant(ldwarf, 10)

# Will need to figure out how to deal with c_1 (supposed to be 2 x angular velocity, i.e., the Coriolis coefficient. Hopefully we don't need c_1 in the
# custom reference framework and will just specify angular_velocity
# If this doesn't work, will need to use override_constants framework

# The "generate transport" scripts will set the transport
# "radial shapes", and the constants c_5, c_6, c_7, c_8, and c_9

the_file = dirname + '/' + fname

print("Writing the atmosphere to %s" %the_file)
print("---------------------------------")
eq.write(the_file)
