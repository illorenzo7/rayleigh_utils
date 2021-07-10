from matplotlib import ticker, colors
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['rapl'])
from common import *
from varprops import *
from plotcommon import *
from rayleigh_diagnostics import GridInfo # for doing averages

# routines to get various Rayleigh slices
def prime(field): # mean along first axis (phi axis)
    shape = np.shape(field)
    shape_collapsed = np.hstack((np.array([1]), shape[1:]))
    return field - np.mean(field, axis=0).reshape(shape_collapsed)

def prime_sph(field, tw): # doesn't work on equatorial slices
    dummy, nt, nr = np.shape(field)
    field_av = np.mean(field, axis=0) # first take the az-avg
    tw_2d = tw.reshape((nt, 1))
    field_av = np.sum(field_av*tw_2d, axis=0)
    return field - field_av.reshape((1, 1, nr))

def get_slice(a, varname, dirname=None, j=0): 
    # gets a basic field associated with Rayleigh object a
    lut = a.lut 
    rr = a.radius
    # first get the appropriate time slice
    vals = a.vals[..., j]
    if vals.ndim == 4 and not hasattr(a, 'lpower'): 
        # Shell_Slice or Meridional_Slice
        cost = a.costheta
        nt = len(cost)
        cost = cost.reshape((1, nt, 1))
    else:
        # note...don't do cylindrical projections for Shell_Spectra!
        # they don't multiply well
        cost = 0.

    # first get root variable name and store any modifiers
    varname, deriv, primevar, sphvar = get_varprops(varname)

    # get sine/cotangent from cosine
    sint = np.sin(np.arccos(cost))
    cott = cost/sint
   
    # shape to make geometric fields
    zero = np.zeros(np.array(np.shape(vals[..., 0])))

    # return the basic field based on the variable name
    if varname in var_indices: 
        # this is really only option for Shell_Spectra...
        the_slice = vals[..., lut[var_indices[varname]]]
    elif varname[-1] in ['l', 'z']: # cylindrical variable
        the_slice_r = vals[..., lut[var_indices[varname[:-1] + 'r']]]
        if deriv:
            the_slice_t = vals[..., lut[var_indices[varname[:-1] + 'T']]]
        else:
            the_slice_t = vals[..., lut[var_indices[varname[:-1] + 't']]]
        if varname[-1] == 'l':
            the_slice = sint*the_slice_r + cost*the_slice_t
        elif varname[-1] == 'z':
            the_slice = cost*the_slice_r - sint*the_slice_t
    elif is_an_int(varname):
        the_slice = vals[..., lut[int(varname)]]
    if primevar:
        the_slice = prime(the_slice)
    elif sphvar:
        gi = GridInfo(dirname + '/grid_info')
        tw = gi.tweights
        the_slice = prime_sph(the_slice, tw)
    del vals # free up memory
    return the_slice

# mollweide + ortho transforms
def mollweide_transform(costheta, clon=0., shrinkage=1., precision=1.e-3): 
    # compute the spherical coordinates
    ntheta = len(costheta)
    nphi = 2*ntheta
    tt = np.arccos(costheta)
    lat = np.pi/2. - tt # these "latitudes" are in radians...
    lon = np.linspace(-np.pi, np.pi, nphi, endpoint=False)

    # compute the beta-angle for the projection 
    # (related to latitude by transcendental equation, which we need
    # to solve iteratively
    beta = np.zeros_like(lat)
    new_beta = np.copy(lat)
    while(np.max(np.abs(new_beta - beta)) > precision):
        beta = new_beta
        new_beta = beta - (2*beta + np.sin(2.*beta) -\
                np.pi*np.sin(lat))/(2.+2.*np.cos(2.*beta))
    # get a "meshgrid" from 1D arrays
    lon, beta = np.meshgrid(lon, beta, indexing='ij')
    xs = 2./np.pi*shrinkage*lon*np.cos(beta)
    ys = shrinkage*np.sin(beta)
    return xs, ys

def ortho_transform(r,lat,lon,lat0=0,lon0=0):
    xs = r*np.cos(lat)*np.sin(lon-lon0)
    ys = r*(np.cos(lat0)*np.sin(lat)-np.sin(lat0)*np.cos(lat)*np.cos(lon-lon0))
    cosc = np.sin(lat0)*np.sin(lat)+np.cos(lat0)*np.cos(lat)*np.cos(lon-lon0) #cosine of angular distance from center of view
    idx = np.where(cosc>=0) #these indices are on front of the globe, the rest should be clipped.
    return xs,ys,idx

# Mollweide plotting routine
plot_moll_kwargs_default = dict({'clon': 0., 'plotlonlines': True, 'lonvals': None, 'plotlatlines': True, 'latvals': np.arange(-60., 90., 30.), 'linewidth': default_lw, 'plotboundary': True})
# change default plotcontours --> False in my_contourf
my_contourf_kwargs_default['plotcontours'] = False
plot_moll_kwargs_default.update(my_contourf_kwargs_default)

def plot_moll(field_orig, costheta, fig, ax, **kwargs):
    kw = update_dict(plot_moll_kwargs_default, kwargs)
    find_bad_keys(plot_moll_kwargs_default, kwargs, 'plot_moll')
    kw_my_contourf = update_dict(my_contourf_kwargs_default, kwargs)
        
    # Shouldn't have to do this but Python is stupid with arrays
    field = np.copy(field_orig)    

    # Get the Mollweide projection coordinates associated with costheta
    xx, yy = mollweide_transform(costheta)

    # shift the field so that the clon is in the ~center of the array
    difflon = 180. - kw.clon # basically difflon is the amount the clon
    # must be shifted to arrive at 180, which is near the center of array
    nphi = 2*len(costheta)
    iphi_shift = int(difflon/360.*nphi)
    field = np.roll(field, iphi_shift, axis=0)

    # make the Mollweide plot
    my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)

    # Draw parallels and meridians, evenly spaced by 30 degrees
    # need some derivative grid info
    tt = np.arccos(costheta)
    lat = np.pi/2. - tt # these "latitudes" are in radians...
    lon = np.linspace(-np.pi, np.pi, 2*len(tt), endpoint=False)
    if kw.lonvals is None:
        kw.lonvals = np.arange(0., 360., 30.)
    
    npoints = 100
    if kw.plotlonlines:
        for lonval in kw.lonvals:
            if lonval == 0.:
                linewidth = 2*kw.linewidth
            else:
                linewidth = kw.linewidth
            # Make sure the plotted meridians are with respect to the clon
            # keep everything in the -180, 180 range
            lon_loc = lonval - kw.clon
            if lon_loc > 180.:
                lon_loc -= 360.
            elif lon_loc < -180.:
                lon_loc += 360.
            lon_loc *= (np.pi/180.)
            imer = np.argmin(np.abs(lon - lon_loc))
            ax.plot(xx[imer, :], yy[imer, :], 'k', linewidth=linewidth)
    if kw.plotlatlines:
        for latval in kw.latvals:
            if latval == 0.: 
                linewidth = 2*kw.linewidth
            else:
                linewidth = kw.linewidth
            ilat = np.argmin(np.abs(lat - latval*np.pi/180.))
            ax.plot(xx[:, ilat], yy[:, ilat], 'k', linewidth=linewidth)

    if kw.plotboundary:
        # Plot outer boundary
        psivals = np.linspace(0, 2*np.pi, 100)
        xvals, yvals = 2.*np.cos(psivals), np.sin(psivals)
        ax.plot(xvals, yvals, 'k', linewidth=1.5*kw.linewidth)

# routine for (l, m) 2D spectra
plot_spec_lm_kwargs_default = dict({'lvals': None, 'mvals': None, 'linewidth': default_lw, 'minmax': None, 'lminmax': None, 'mminmax': None,\
    # more cbar stuff
    'plotcbar': True, 'cmap': None, 'norm': None, 'linear': False, 'units': '', 'fontsize': default_labelsize})
plot_spec_lm_kwargs_default.update(add_cbar_kwargs_default)

def plot_spec_lm(field, fig, ax, **kwargs):
    kw = update_dict(my_contourf_kwargs_default, kwargs)
    kw_add_cbar = update_dict(add_cbar_kwargs_default, kwargs)
    find_bad_keys(plot_spec_lm_kwargs_default, kwargs, 'plot_spec_lm')

    # make sure Python does not modify any of the arrays it was passed
    field = np.copy(field)

    # full (l, m) grid:
    nell, nm = np.shape(field)
    lvals_all = np.arange(nell)
    mvals_all = np.arange(nm)

    if not kw.lminmax is None:
        il1 = np.argmin(np.abs(lvals_all - kw.lminmax[0]))
        il2 = np.argmin(np.abs(lvals_all - kw.lminmax[1]))
    else:
        il1, il2 = 0, nell - 1

    if not kw.mminmax is None:
        im1 = np.argmin(np.abs(mvals_all - kw.mminmax[0]))
        im2 = np.argmin(np.abs(mvals_all - kw.mminmax[1]))
    else:
        im1, im2 = 0, nm - 1

    # now adjust everything by the (l, m) range we want
    lvals = lvals_all[il1:il2+1]
    mvals = mvals_all[im1:im2+1]
    field = field[il1:il2+1, im1:im2+1]

    lvals_2d, mvals_2d = np.meshgrid(lvals, mvals, indexing='ij')
    lvals_2d, mvals_2d = xy_grid(lvals_2d, mvals_2d)

    # Get minmax, if not specified
    if kw.minmax is None:
        field_not0 = np.copy(field)
        # power gets wierd (close to 0?) at the two
        # most extreme l-values
        if il1 == 0: 
            field_not0 = field_not0[1:, :]
        if il2 == nell - 1: 
            field_not0 = field_not0[:-1, :]
        field_not0 = field_not0[field_not0 != 0.]
        if kw.linear: # NOT the default...
            kw.minmax = kw_add_cbar.minmax =\
                contourf_minmax(field_not0, posdef=True)
            kw_add_cbar.posdef = True
        else:
            kw.minmax = kw_add_cbar.minmax =\
                contourf_minmax(field_not0, logscale=True, buff_ignore1=None, buff_ignore2=None)
            kw_add_cbar.logscale = True
   
    if kw.norm is None and not kw.linear: # the default
        kw.norm = colors.LogNorm(vmin=kw.minmax[0], vmax=kw.minmax[1])
    
    if kw.cmap is None:
        kw.cmap = 'jet'
    im = plt.pcolormesh(lvals_2d, mvals_2d, field, cmap=kw.cmap, norm=kw.norm)  

    # now deal with color bar, if one is desired
    if kw.plotcbar:
        add_cbar(fig, ax, im, **kw_add_cbar)
