import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import Shell_Avgs, GridInfo
from common import *

def get_numbers_input(dirname, r1='rmin', r2='rmax'):
    di = dotdict()
    rotation = get_parameter(dirname, 'rotation')
    magnetism = get_parameter(dirname, 'magnetism')

    # get non-rotating, non-magnetic numbers first:

    # aspect ratio
    r1, r2 = interpret_rvals(dirname, np.array([r1, r2]))
    di.r1, di.r2 = r1, r2
    di.aspect = r2/r1

    # density contrast
    eq = get_eq(dirname)
    gi = get_grid_info(dirname)
    rr = gi.rr
    ir1 = np.argmin(np.abs(rr - r1))
    ir2 = np.argmin(np.abs(rr - r2))
    di.dc = eq.rho[ir1]/eq.rho[ir2]

    # Prandtl number


    return di
