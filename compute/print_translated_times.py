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
val_iter = clas['val_iter']
di = translate_times(val_iter, dirname)

val_iter = di['val_iter']
val_unit = val_tdt = di['val_tdt']
label = 'TDTs'
val_sec = di['val_sec']
if rotation:
    val_unit = di.val_prot
    label = 'rotations'

print ('%08i iter =\n%.3f %s' %(val_iter, val_unit, label))
