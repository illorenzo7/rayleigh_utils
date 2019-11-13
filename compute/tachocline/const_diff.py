# Author: Loren Matilsky
# Created: 11/12/2019
#
# Purpose: modify a binary file (default name custom_reference_binary) 
# to contain diffusive profiles that are constant
#
# Reference state (with grid info) must already be present in
# custom_reference_binary
#
# Parameters: output_dir (first argument), 
# command line options:
#
# -nutop
# default 3 x 10^12 c.g.s
#
# -kappatop
# default 3 x 10^12 c.g.s
#
# -etatop
# default 3 x 10^12 c.g.s
#
# -mag
# Whether magnetism is True or False, default False (hydro)

import numpy as np
import sys, os

sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])

from reference_tools import equation_coefficients

# Set default constants
nutop = 3.0e12
kappatop = 3.0e12
etatop = 3.0e12
mag = False

# Get directory to save binary files for reference state and heating
dirname = sys.argv[1]
fname = 'custom_reference_binary'

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-nutop':
        nutop = float(args[i+1])
    elif arg == '-kappatop':
        kappatop = float(args[i+1])
    elif arg == '-etatop':
        etatop = float(args[i+1])
    elif arg == '-fname':
        fname = args[i+1]
    elif arg == '-mag':
        mag = True

# If hydro, better make sure whatever multiplies eta in energy/induction
# is zero... but "had we better?" does eta != 0 cause crashes?
if not mag:
    etatop = 0.0

# Open and read the hopefully already existing reference file!
eq = equation_coefficients()
the_file = dirname + '/' + fname
eq.read(the_file)
nr = eq.nr

# Now write to file using the equation_coefficients framework
print("Setting f_3, f_5, f_7, f_11, f_12, and f_13")
eq.set_function(nutop*np.ones(nr), 3)
eq.set_function(kappatop*np.ones(nr), 5)
eq.set_function(etatop*np.ones(nr), 7)

eq.set_function(np.zeros(nr, dtype='float'), 11)
eq.set_function(np.zeros(nr, dtype='float'), 12)
eq.set_function(np.zeros(nr, dtype='float'), 13)

print("Setting c_5, c_6, c_7, c_8, and c_9")
if not mag:
    print ("magnetism = False, so c_7 and c_9 = 0")

eq.set_constant(1.0, 5) # multiplies viscous force
eq.set_constant(1.0, 6) # multiplies thermal diffusion term
eq.set_constant(1.0, 7) # multiplies eta in induction equation
eq.set_constant(1.0, 8) # multiplies viscous heating term
eq.set_constant(1.0/4.0/np.pi, 9) # multiplies magnetic diffusion term

# Will need to figure out how to deal with c_1 (supposed to be 2 x angular velocity, i.e., the Coriolis coefficient. Hopefully we don't need c_1 in the
# custom reference framework and will just specify angular_velocity
# If this doesn't work, will need to use override_constants framework

# c_10 will be set in the "generate heating" scripts

print("Writing the atmosphere to %s" %the_file)
print("---------------------------------")
eq.write(the_file)
