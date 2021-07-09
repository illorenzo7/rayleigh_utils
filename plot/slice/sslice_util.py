from matplotlib import ticker, colors
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['rapl'])
from common import *
from plotcommon import *

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
