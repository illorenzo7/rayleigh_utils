# Author: Loren Matilsky
# Created: 12/10/2019
#
# Purpose: generate a binary file (for Rayleigh to read) that contains
# a reference state, using the old "reference" and "transport" files 

# Parameters: 
# input_dir (first argument; where the reference/transport files are)
# output_dir (second argument; where the custom_reference_binary file
    # WILL go

import numpy as np
import sys, os
import struct
from scipy.integrate import cumtrapz

sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])

from reference_tools import equation_coefficients
from rayleigh_diagnostics import ReferenceState, TransportCoeffs

cp = 3.5e8

# Get directory to save binary files for reference state and heating
inputdir = sys.argv[1]
outputdir = sys.argv[2]

# Name for output custom reference state
fname = 'custom_reference_binary'

ref = ReferenceState(inputdir + '/reference', '')
trans = TransportCoeffs(inputdir + '/transport', '')

print ("Read reference stuff from " + inputdir + '/reference')
print ("and " + inputdir + '/transport')

# Radial grid
nr = ref.nr
rr = ref.radius

try: 
    eta = trans.eta # if this succeeds, magnetism must be "on"
    dlneta = trans.dlneta
    mag = True
    print ("Found eta(r) and dlneta(r) in transport, setting mag = True")
except:
    eta = np.zeros(nr)
    dlneta = np.zeros(nr)
    mag = False
    print("No eta(r) or dlneta(r) in transport, setting mag = False")

# Calculate the heating function f_6 and luminosity c_10
rhot = ref.density*ref.temperature
Q = ref.heating*rhot

lum = -4.*np.pi*cumtrapz(rr**2*Q, rr)[-1]
print('luminosity from ref.heating using cumtrapz is %1.3e' %lum)

# Set the heating coefficients
f6 = Q/lum
c10 = lum

# Now write to file using the equation_coefficients framework
eq = equation_coefficients(rr)
# Only set c_4 = 1/(4*pi) if mag = True

print("Setting all functions")
eq.set_function(ref.density, 1)
buoy = ref.density*ref.gravity/cp
eq.set_function(buoy, 2)
eq.set_function(trans.nu, 3)
eq.set_function(ref.temperature, 4)
eq.set_function(trans.kappa, 5)
eq.set_function(f6, 6)
eq.set_function(eta, 7)
eq.set_function(ref.dlnrho, 8)
eq.set_function(ref.d2lnrho, 9)
eq.set_function(ref.dlnt, 10)
eq.set_function(trans.dlnu, 11)
eq.set_function(trans.dlnkappa, 12)
eq.set_function(dlneta, 13)
eq.set_function(ref.dsdr, 14)

print("Setting all constants except c_1")
print("c_1 should be overridden in main_input for rotating models")
eq.set_constant(1.0, 2) # multiplies buoyancy

eq.set_constant(1.0, 3) # multiplies pressure grad.

if mag:
    print("magnetism = True, so setting c_4 = 1/(4*pi)")
    eq.set_constant(1.0/4.0/np.pi, 4)
else:
    print("magnetism = False, so setting c_4 = 0.0")
    print("Should also set c_7 and c_9 to 0, but didn't for other ")
    print("custom reference state stuff")
    eq.set_constant(0.0, 4)

eq.set_constant(1.0, 5) # multiplies viscous force

eq.set_constant(1.0, 6) # multiplies thermal dissipation

eq.set_constant(1.0, 7) # multiplies nu(r) in induction eqn.
    # SHOULD TECHNICALLY BE ZERO IF MAG = FALSE

eq.set_constant(1.0, 8) # multiplies viscous dissipation

eq.set_constant(1.0/4.0/np.pi, 9) # multiplies Ohmic dissipation
    # SHOULD TECHNICALLY BE ZERO IF MAG = FALSE

eq.set_constant(c10, 10) # luminosity

# Now write equation coefficients to file
the_file = outputdir + '/' + fname

print("Writing the atmosphere to %s" %the_file)
print("---------------------------------")
eq.write(the_file)
