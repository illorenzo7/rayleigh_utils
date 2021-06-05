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
di_trans = clas['di_trans']

val_iter = di_trans['val_iter']
val_unit = di_trans['val_unit']
simple_label = di_trans['simple_label']
val_sec = di_trans['val_sec']

print ('%08i iter =\n%.3f %s' %(val_iter, val_unit, simple_label))
