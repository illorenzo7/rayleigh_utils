from matplotlib import ticker, colors
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['raco'] + '/quantities_util')
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['rapl'])
from common import *
from varprops import *
from plotcommon import *
from rayleigh_diagnostics import GridInfo # for doing averages

# default fig dimensions
generic_slice_fig_dimensions = dict({'sub_height_inches': 3., 'sub_margin_left_inches': 1/4, 'sub_margin_right_inches': 1/4, 'sub_margin_top_inches': 1, 'sub_margin_bottom_inches': 1/2, 'margin_top_inches': default_margin})

moll_fig_dimensions = dict({'sub_aspect': 1/2})
ortho_fig_dimensions = dict({'sub_aspect': 1})

moll_fig_dimensions.update(generic_slice_fig_dimensions)
ortho_fig_dimensions.update(generic_slice_fig_dimensions)

eq_fig_dimensions = ortho_fig_dimensions

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

    # first test if varname is valid and if it's a basic variable
    basic = is_basic(varname)

    # first get the appropriate time slice
    vals = a.vals[..., j]
    lut = a.lut 
    rr = a.radius

    if basic:
        # gets a basic field associated with Rayleigh object a
        if vals.ndim == 4 and not hasattr(a, 'lpower'): 
            # Shell_Slice or Meridional_Slice
            cost = a.costheta
            nt = len(cost)
            cost = cost.reshape((1, nt, 1))
        else: #an Equatorial_Slices or Shell_Spectra
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

        # for abs om and pv, right now this is just the vorticity
        if 'absom' in varname or 'pv' in varname: 
            # get the planetary vorticity
            eq = get_eq(dirname)
            if varname[-1] == 'z': # easiest
                the_slice += 2*eq.om0
            elif varname[-1] == 'r': 
                the_slice += 2*eq.om0*cost
            elif varname[-1] == 't':
                the_slice -= 2*eq.om0*sint
            # for lambda and phi, there is no planetary vorticity 
            # contribution

        if 'pv' in varname: # need to divide by H
            the_slice /= compute_axial_H(rr, sint)

        if primevar:
            the_slice = prime(the_slice)
        elif sphvar:
            gi = GridInfo(dirname + '/grid_info')
            tw = gi.tweights
            the_slice = prime_sph(the_slice, tw)
        del vals # free up memory

        return the_slice
    else:
        if '+' in varname or '=' in varname:
            signature = [1]
            for char in varname:
                if char == '+':
                    signature.append(1)
                elif char == '=':
                    signature.append(-1)

            the_slice = np.zeros_like(vals[..., 0])
            count = 0
            for subvar in varname.split('+'):
                for subsubvar in subvar.split('='):
                    the_slice += get_slice(a, subsubvar, dirname, j)*signature[count]
                    count += 1
        elif '*' in varname or '/' in varname:
            signature = [1]
            for char in varname:
                if char == '*':
                    signature.append(1)
                elif char == '/':
                    signature.append(-1)

            the_slice = np.ones_like(vals[..., 0])
            count = 0
            for subvar in varname.split('*'):
                for subsubvar in subvar.split('/'):
                    the_slice *= (get_slice(a, subsubvar, dirname, j))**signature[count]
                    count += 1
        return the_slice

# mollweide + ortho transforms
def mollweide_transform(costheta, shrinkage=1., precision=1.e-3): 
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
    xx = 2./np.pi*shrinkage*lon*np.cos(beta)
    yy = shrinkage*np.sin(beta)
    return xx, yy

def ortho_transform(costheta, clat=20., shrinkage=1.):
    # compute the spherical coordinates
    ntheta = len(costheta)
    nphi = 2*ntheta
    tt = np.arccos(costheta)
    lat = np.pi/2. - tt # these "latitudes" are in radians...
    lon = np.linspace(-np.pi, np.pi, nphi, endpoint=False)
    clat = np.pi/180*clat # convert clat, degrees --> radians

    # get a "meshgrid" from 1D arrays
    lon, lat = np.meshgrid(lon, lat, indexing='ij')
   
    # do ortho projection
    xx = shrinkage*np.cos(lat)*np.sin(lon)
    yy = shrinkage*(np.cos(clat)*np.sin(lat) -\
        np.sin(clat)*np.cos(lat)*np.cos(lon))
    cosc = np.sin(clat)*np.sin(lat)+np.cos(clat)*np.cos(lat)*np.cos(lon) #cosine of angular distance from center of view
    igood = np.where(cosc>=0) #these indices are on front of the globe, the rest should be clipped.
    ibad = np.where(cosc<0)
    return xx,yy,igood,ibad

# Mollweide ortho ortho plotting routine
plot_moll_or_ortho_kwargs_default = dict({'clon': 0., 'clat': 20., 'shrinkage': 1., 'plotlonlines': True, 'lonvals': np.arange(0., 360., 30.), 'plotlatlines': True, 'latvals': np.arange(-60., 90., 30.), 'linewidth': default_lw, 'plotboundary': True, 'ortho': False})
plot_moll_or_ortho_kwargs_default.update(my_contourf_kwargs_default)
plot_moll_or_ortho_kwargs_default['plotcontours'] = False

def plot_moll_or_ortho(field, costheta, fig, ax, **kwargs):
    kw = update_dict(plot_moll_or_ortho_kwargs_default, kwargs)
    find_bad_keys(plot_moll_or_ortho_kwargs_default, kwargs, 'plot_moll_or_otho')
    # change default plotcontours --> False in my_contourf
    tmp = my_contourf_kwargs_default.copy()
    tmp['plotcontours'] = False
    kw_my_contourf = update_dict(tmp, kwargs)
        
    # Shouldn't have to do this but Python is stupid with arrays
    field = np.copy(field)    

    # shift the field so that the clon is in the ~center of the array
    difflon = 180. - kw.clon # basically difflon is the amount the clon
    # must be shifted to arrive at 180, which is near the center of array
    nphi = 2*len(costheta)
    iphi_shift = int(difflon/360.*nphi)
    field = np.roll(field, iphi_shift, axis=0)

    # Get the projection coordinates associated with costheta
    if kw.ortho:
        xx, yy, igood, ibad = ortho_transform(costheta, clat=kw.clat, shrinkage=kw.shrinkage)

        field[ibad] = np.nan
        xx_masked = np.copy(xx)
        yy_masked = np.copy(yy)
        xx_masked[ibad] = np.nan
        yy_masked[ibad] = np.nan
    else:
        xx, yy = mollweide_transform(costheta)
        xx_masked = np.copy(xx)
        yy_masked = np.copy(yy)

    # make the color plot plot
    my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)

    # Draw parallels and meridians, evenly spaced by 30 degrees
    # need some derivative grid info
    tt = np.arccos(costheta)
    lat = np.pi/2. - tt # these "latitudes" are in radians...
    lon = np.linspace(-np.pi, np.pi, 2*len(tt), endpoint=False)
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
            ax.plot(xx_masked[imer, :], yy_masked[imer, :], 'k', linewidth=linewidth)
    if kw.plotlatlines:
        for latval in kw.latvals:
            if latval == 0.: 
                linewidth = 2*kw.linewidth
            else:
                linewidth = kw.linewidth
            ilat = np.argmin(np.abs(lat - latval*np.pi/180.))
            ax.plot(xx_masked[:, ilat], yy_masked[:, ilat], 'k', linewidth=linewidth)

    if kw.plotboundary:
        # Plot outer boundary
        psivals = np.linspace(0, 2*np.pi, 100)
        xvals, yvals = 2.*np.cos(psivals), np.sin(psivals)
        if kw.ortho:
            xvals /= 2.
        ax.plot(xvals, yvals, 'k', linewidth=1.5*kw.linewidth)

# equatorial slice plotting routine
plot_eq_kwargs_default = dict({'clon': 0., 'plotlonlines': True, 'lonvals': np.arange(0., 360., 60.), 'linewidth': default_lw, 'plotboundary': True})
plot_eq_kwargs_default.update(my_contourf_kwargs_default)
plot_eq_kwargs_default['plotcontours'] = False

def plot_eq(field, rr, fig, ax, **kwargs):
    kw = update_dict(plot_eq_kwargs_default, kwargs)
    find_bad_keys(plot_eq_kwargs_default, kwargs, 'plot_eq')
    # change default plotcontours --> False in my_contourf
    tmp = my_contourf_kwargs_default.copy()
    tmp['plotcontours'] = False
    kw_my_contourf = update_dict(tmp, kwargs)
        
    # Shouldn't have to do this but Python is stupid with arrays
    field = np.copy(field)    

    # shift the field so that the clon is in the ~center of the array
    difflon = 180. - kw.clon # basically difflon is the amount the clon
    # must be shifted to arrive at 180, which is near the center of array
    nphi, nr = np.shape(field)
    iphi_shift = int(difflon/360.*nphi)
    field = np.roll(field, iphi_shift, axis=0)

    # get the projection coordinates
    # extend phi (avoid little wedge) ...
    nphi += 1
    phi = np.linspace(-np.pi, np.pi, nphi)
    lon = 180*phi/np.pi
    cosphi = np.cos(phi)
    sinphi = np.sin(phi)

    # extend the field
    field = np.vstack((field, field[0].reshape((1, nr))))

    # 2 D r and phi
    rr_2d = rr.reshape((1, nr))
    phi_2d = phi.reshape((nphi, 1))

    # cartesian grid # put phi=0 lower center
    rmax = np.max(rr)
    xx = rr_2d*np.sin(phi_2d)/rmax
    yy = -rr_2d*np.cos(phi_2d)/rmax

    # make the color plot
    my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)

    # longitude lines
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
            ilon = np.argmin(np.abs(lon - lon_loc))
            ax.plot(xx[ilon, :], yy[ilon, :], 'k', linewidth=linewidth)

    if kw.plotboundary:
        # plot boundaries
        ax.plot(xx[:, 0], yy[:, 0], 'k', linewidth=1.5*kw.linewidth)
        ax.plot(xx[:, -1], yy[:, -1], 'k', linewidth=1.5*kw.linewidth)
