import sys, os
sys.path.append(os.environ['rapl'])
sys.path.append(os.environ['raco'])
from common import *
from lut import *

def slice_ndarray(arr, axis, index):
    arange = list(np.arange(arr.ndim))
    shapeout = list(np.shape(arr))
    shapeout[axis] = 1
    permute = [arange.pop(axis)] + arange # permute so the array has the 
        # desired axis first
    outarr = np.transpose(arr, axes=permute)[index]
    return outarr.reshape(shapeout)

def derive_quantity(dirname, vals, lut, quantity_or_index, eqslice=False, timeaxis=False):
    # assumes vals of shape
    # (nphi, ntheta, nr, nq)
    # (ntheta, nr, nq)
    # (nr, nq)
    # or 
    # (nphi, nr, nq) (this one is special (eqslice) and must be specified)

    # parse the quantity (could be user (Loren) - defined
    index, quantity = parse_quantity(quantity_or_index)

    # get reference state
    eq = get_eq(dirname)

    shape_vals = list(np.shape(vals))
    shape_needed = np.ones(vals.ndim, 'int') # shape of ref quantity should be 
        # one less than vals shape after grabbing a specific quantity
    
    if timeaxis: # last dims are (..., nr, nq, niter)
        radaxis = vals.ndim - 3
        qaxis = vals.ndim - 2
    else: # last dims are (..., nr, nq)
        radaxis = vals.ndim - 2
        qaxis = vals.ndim - 1
    shape_needed[radaxis] = shape_vals[radaxis]

    rho = eq.density.reshape(shape_needed)
    T = eq.temperature.reshape(shape_needed)
    dSdr = eq.dsdr.reshape(shape_needed)
    
    if quantity == 'advref':
        vr = slice_ndarray(vals, qaxis, lut[1])
        return rho*T*dSdr*vr   
    if quantity == 'rhov_r':
        vr = slice_ndarray(vals, qaxis, lut[1])
        return rho*vr
    if quantity == 'rhov_theta':
        vt = slice_ndarray(vals, qaxis, lut[2])
        return rho*vt
    if quantity == 'rhov_phi':
        vp = slice_ndarray(vals, qaxis, lut[3])
        return rho*vp
