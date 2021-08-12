import numpy as np
import scipy.special as ss
from compute_grid_info import compute_theta_grid

def to_spec(field, nell=None, nm=None, phi=None, tt=None, tw=None, pw=None):
    nphi, nt = np.shape(field)
    # check nphi = 2*nt
    if nphi != 2*nt:
        print ("nphi != 2*nt")
        print ("Houston, we may have a problem!")

    # dealias by default
    if nell is None:
        nell = 2*nt//3
    if nm is None:
        nm = 2*nphi//3

    # make the grid and weights if not provided
    tt_default, tw_default = compute_theta_grid(nt)
    if tw is None:
        tw = tw_default
    tw = tw.reshape((1, nt))

    if pw is None:
        pw = np.zeros(nphi) + 1/nphi
    pw = pw.reshape((nphi, 1))

    if tt is None:
        tt = tt_default

    if phi is None:
        phi = np.linspace(0, 2*np.pi, nphi, endpoint=False)

    phi_2d, tt_2d = np.meshgrid(phi, tt, indexing='ij')
    phiflat, ttflat = phi_2d.flatten(), tt_2d.flatten()

    # now make the spectral array
    field_spec = np.zeros((nell, nm), dtype='complex')
    for im in range(nm):
        for il in range(im, nell):
            ylmstar = ss.sph_harm(-im, il, phiflat, ttflat).reshape((nphi, nt))
            field_spec[il, im] = np.sum(np.sum(field*ylmstar*pw*tw, axis=0), axis=0)/np.sqrt(4*np.pi)
        if im > 0:
            field_spec[:, im] /= np.sqrt(2)

    return field_spec


