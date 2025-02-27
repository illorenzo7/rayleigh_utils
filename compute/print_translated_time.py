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

# only need one cla here (iter, t_omega, or t_kappa)
for key, val in clas.items():
    tvals = make_array(val) # make this an array
    ntvals = len(tvals)
    translate_from = key

    val_iter = np.zeros(ntvals, dtype='int')
    val_simt = np.zeros(ntvals, dtype='float')
    val_tkappa = np.zeros(ntvals, dtype='float')
    if clas0.rotation:
        val_tomega = np.zeros(ntvals, dtype='float')

    for i in range(ntvals):
        tval = tvals[i] # we have multiple tvals
        di = translate_times(tval, dirname, translate_from=translate_from)
        val_iter[i] = di.val_iter
        val_simt[i] = di.val_simt
        val_tkappa[i] = di.val_tkappa
        if clas0.rotation:
            val_tomega[i] = di.val_tomega

    print ("val_iter =", arr_to_str(val_iter, '%08i', nobra=True))
    print ("val_simt =", arr_to_str(val_simt, '%.2f', nobra=True))
    print ("val_tkappa =", arr_to_str(val_tkappa, '%.3f', nobra=True))
    if clas0.rotation:
        print ("val_tomega =", arr_to_str(val_tomega, '%.2f', nobra=True))
