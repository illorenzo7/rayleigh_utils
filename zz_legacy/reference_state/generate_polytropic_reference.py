# Loren Matilsky
# 06/27/2020
# computes polytropic reference analytically and writes to 
# "custom_reference binary"
# need --ro (location of outside radius, as fraction of solar radius
# everything else is default:
# --ri 0.719
# --poly_n 1.5
# --nrho 3      OR
# --r1 0.9467 (ro = 6.586d10)

import numpy as np
import sys, os
from arbitrary_atmosphere import arbitrary_atmosphere

sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
#sys.path.append(os.environ['raco'] + '/reference')

from reference_tools import equation_coefficients
from common import *
from cla_util import *
from polytrope import *

# Get directory to save binary files for reference
clas0, clas = read_clas_raw(sys.argv)
dirname = clas0['dirname']
fname = 'custom_reference_binary'

# overwrite defaults
kw = update_dict(compute_polytrope2_kwargs_default, clas)

# check for bad keys
find_bad_keys(compute_polytrope2_kwargs_default, clas, clas0['routinename'], justwarn=True)

print(buff_line)
print ("computing polytropic reference:")
print_dict(kw)
di_ref = compute_polytrope2(**clas)
rho = di_ref['rho']
dlnrho = di_ref['dlnrho']
d2lnrho = di_ref['d2lnrho']
T = di_ref['T']
dlnT = di_ref['dlnT']
dSdr = di_ref['dSdr']
rr = di_ref['rr']
g = G*kw.mstar/rr**2

# Now write to file using the equation_coefficients framework
print(buff_line)
the_file = dirname + '/' + fname
print("Writing the atmosphere to %s" %the_file)
eq = equation_coefficients(rr)

# Set only the thermodynamic functions/constants in this routine
# In other routines, we can set the heating and transport coefficients
print("Setting f_1, f_2, f_4, f_8, f_9, f_10, and f_14")
eq.set_function(rho, 1)
buoy = rho*g/c_P
eq.set_function(buoy, 2)
eq.set_function(T, 4)
eq.set_function(dlnrho, 8)
eq.set_function(d2lnrho, 9)
eq.set_function(dlnT, 10)
eq.set_function(dSdr, 14)

print("Setting c_2, c_3, c_7, and c_8")
eq.set_constant(1.0, 2) # multiplies buoyancy
eq.set_constant(1.0, 3) # multiplies pressure grad.
eq.set_constant(1.0, 8) # multiplies viscous heating

eq.write(the_file)
print (buff_line)
