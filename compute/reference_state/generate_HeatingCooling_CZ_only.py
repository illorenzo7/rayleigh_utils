# Author: Loren Matilsky
# Created: 11/29/2023
#
# Purpose: generate a binary file (for Rayleigh to read) that contains
# a heating profile Q(r) = c_10*f_6
# that enforces a convectively unstable layer 
#
# parameters:
# beta: aspect ratio of CZ
# deltain: width of bottom heating layer
# deltac: width of top cooling layer

import numpy as np
import sys, os

sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['raco'] + '/reference_state')

from reference_tools import equation_coefficients
from common import *
from cla_util import *
from Jupiter_Heating_CZ_WL import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas_raw(args)
dirname = clas0['dirname']

# Set default kwargs
kw_default = dotdict(dict({'nr': nr_default, 'beta': beta_default, 'deltah': deltah_default, 'deltac': deltac_default}))

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# overwrite defaults
kw = update_dict(kw_default, clas)

# compute heating function
rin, rout, rr, heating = compute_heating_CZ_only(**kw)

# print what we computed
print(buff_line)
print("Computed heating/cooling function for CZ")
print("nr         : %i" %kw.nr) 
print("rin        : %.3f" %rin)
print("rout       : %.3f" %rout)
print(buff_line)
print("beta       : %.3f" %kw.beta) 
print("deltah     : %.3f" %kw.deltah)
print("deltac     : %.3f" %kw.deltac)
print(buff_line)

# Now write to file using the equation_coefficients framework
eq = equation_coefficients(rr)

# Set only c_10*f_6
print("Setting c_10 and f_6")
eq.set_function(heating, 6)
c10 = 4.*np.pi*kw.beta**2/(1.-kw.beta)**2
eq.set_constant(c10, 10) 
print("c_10       : %1.3e" %c10)
print(buff_line)

fname = "HeatingCoolingCZ_beta%.3f_deltah%.3f_deltac%.3f" %(kw.beta, kw.deltah, kw.deltac)

the_file = dirname + '/' + fname

print("Writing the atmosphere to %s" %the_file)
print(buff_line)
eq.write(the_file)
