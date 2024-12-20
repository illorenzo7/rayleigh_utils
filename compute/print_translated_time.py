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

# only need one cla here (iter, prot, or tdt)
for key, val in clas.items():
    tvals = make_array(val) # make this an array
    ntvals = len(tvals)
    translate_from = key

    val_iter = np.zeros(ntvals, dtype='int')
    val_simt = np.zeros(ntvals, dtype='float')
    val_tdt = np.zeros(ntvals, dtype='float')
    if clas0.rotation:
        val_prot = np.zeros(ntvals, dtype='float')

    for i in range(ntvals):
        tval = tvals[i] # we have multiple tvals
        di = translate_times(tval, dirname, translate_from=translate_from)
        val_iter[i] = di.val_iter
        val_simt[i] = di.val_simt
        val_tdt[i] = di.val_tdt
        if clas0.rotation:
            val_prot[i] = di.val_prot

    print ("val_iter =", arr_to_str(val_iter, '%08i', nobra=True))
    print ("val_simt =", arr_to_str(val_simt, '%.2f', nobra=True))
    print ("val_tdt =", arr_to_str(val_tdt, '%.3f', nobra=True))
    if clas0.rotation:
        print ("val_prot =", arr_to_str(val_prot, '%.2f', nobra=True))
