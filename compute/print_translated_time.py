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

    iters = np.zeros(ntvals, dtype='int')
    times = np.zeros(ntvals, dtype='float')
    tkappa = np.zeros(ntvals, dtype='float')
    if clas0.rotation:
        tomega = np.zeros(ntvals, dtype='float')

    for i in range(ntvals):
        tval = tvals[i] # we have multiple tvals
        di = translate_times(tval, dirname, translate_from=translate_from)
        iters[i] = di.iter
        times[i] = di.time
        tkappa[i] = di.tkappa
        if clas0.rotation:
            tomega[i] = di.tomega

    print ("iters =", arr_to_str(iters, '%08i', nobra=True))
    print ("times =", arr_to_str(times, '%1.4e', nobra=True))
    print ("tkappa =", arr_to_str(tkappa, '%1.4e', nobra=True))
    if clas0.rotation:
        print ("tomega =", arr_to_str(tomega, '%1.4e', nobra=True))
