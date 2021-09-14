import sys, os
sys.path.append(os.environ['rapl'])
sys.path.append(os.environ['raco'])
from common import *
from lut import *

def derived_quantity(dirname, vals, lut, quantity_or_index, eqslice=False):
    # assumes vals of shape
    # (nphi, ntheta, nr, nq)
    # (ntheta, nr, nq)
    # (nr, nq)
    # or 
    # (nphi, nr, nq) (this one is special (eqslice) and must be specified)

    index, quantity = parse_quantity(quantity_or_index)

    # some quantities will be listed as strings to be derived
    eq = get_eq(dirname)

    shape_vals = np.shape(vals)
    shape_needed = np.ones(vals.ndim - 1, 'int')
    shape_needed[-1] = shape_vals[-2]

    rho = eq.density.reshape(shape_needed)
    T = eq.temperature.reshape(shape_needed)
    dSdr = eq.dsdr.reshape(shape_needed)

    if quantity == 'advref':
        return rho*T*dSdr*vals[:, :, lut[1]]
    if quantity == 'rhov_r':
        return rho*vals[:, :, lut[1]]
    if quantity == 'rhov_theta':
        return rho*vals[:, :, lut[2]]
    if quantity == 'rhov_phi':
        return rho*vals[:, :, lut[3]]
