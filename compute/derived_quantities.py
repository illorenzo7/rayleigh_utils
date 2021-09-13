import sys, os
sys.path.append(os.environ['rapl'])
sys.path.append(os.environ['raco'])
from common import *

def derived_azav(dirname, vals, lut, qval):
    # some quantities will be listed as strings to be derived
    eq = get_eq(dirname)
    nt, nr, dummy = np.shape(vals)
    rho = eq.density.reshape((1, nr))
    T = eq.temperature.reshape((1, nr))
    dSdr = eq.dsdr.reshape((1, nr))
    if qval == 'advref':
        return rho*T*dSdr*vals[:, :, lut[1]]
