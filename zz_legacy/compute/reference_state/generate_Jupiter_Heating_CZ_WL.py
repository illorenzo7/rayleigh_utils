# Author: Loren Matilsky
# Created: 11/28/2023
#
# Purpose: generate a binary file (for Rayleigh to read) that contains
# a heating profile Q(r) = c_10*f_6
# that enforces a convectively unstable layer below a (slightly) 
# stable layer
#
# parameters:
# alpha: aspect ratio of WL-to-CZ
# beta: aspect ratio of full layer
# delta1: width of bottom heating layer
# delta2: width of top cooling layer
# fluxratio: ratio of stable flux to unstable flux

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
kw_default = dotdict(dict({'nr': nr_default,  'alpha': alpha_default, 'beta': beta_default, 'delta1': delta1_default, 'delta2': delta2_default, 'fluxratio': fluxratio_default}))

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# overwrite defaults
kw = update_dict(kw_default, clas)

# compute heating function
rin, rout, r0, rr, heating = generate_heating_CZ_WL(**kw)

# print what we computed
print(buff_line)
print("Computed heating function for WL atop CZ")
print("nr         : %i" %kw.nr) 
print("rin        : %.3f" %rin)
print("r0         : %.3f" %r0)
print("rout       : %.3f" %rout)
print(buff_line)
print("alpha      : %1.3e" %kw.alpha) 
print("beta       : %1.3e" %kw.beta) 
print("delta1     : %1.3e" %kw.delta1)
print("delta2     : %1.3e" %kw.delta1)
print("fluxratio  : %1.3e" %kw.fluxratio)
print(buff_line)

# Now write to file using the equation_coefficients framework
eq = equation_coefficients(rr)

# Set only c_10*f_6
print("Setting c_10 and f_6")
eq.set_function(heating, 6)
eq.set_constant(1.0, 10) 

fname = "customfile_alpha%1.3e_beta%1.3e_delta1%1.3e_delta2%1.3e_fluxratio%1.3e" %(kw.alpha, kw.beta, kw.delta1, kw.delta2, kw.fluxratio)

the_file = dirname + '/' + kw.fname

print("Writing the atmosphere to %s" %the_file)
print("---------------------------------")
eq.write(the_file)
