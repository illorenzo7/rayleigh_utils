# Author: Loren Matilsky
# Date created: well before 05/06/2019
# Prints various timescales associated with a given Rayleigh simulation
# directory, including how long it was run, the thermal diffusion time, 
# rotation period, etc.

import sys
import numpy as np
from common import *
from cla_util import *

# Get CLAs
args = sys.argv 
clas0, clas = read_clas(args)
dirname = clas0['dirname']
rotation = clas0.rotation

rotation = clas0.rotation

# only need one clas here (iter, prot, tdt, sec)
for key, val in clas.items():
    translate_from = key
    tval = val

di = translate_times(tval, dirname, translate_from=translate_from)
for key, val in di.items():
    print ("%8s = %1.3e" %(key,val))
