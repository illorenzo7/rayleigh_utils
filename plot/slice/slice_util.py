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
generic_slice_fig_dimensions = dict({'sub_height_inches': 3., 'sub_margin_left_inches': default_margin, 'sub_margin_right_inches': default_margin, 'margin_left_inches': 0., 'margin_right_inches': 0.,'sub_margin_top_inches': 1, 'sub_margin_bottom_inches': 3/8, 'margin_top_inches': default_margin})

moll_fig_dimensions = dict({'sub_aspect': 1/2})
molllonav_fig_dimensions = dict({'sub_aspect': 1/2, 'sub_margin_right_inches': 2.})
ortho_fig_dimensions = dict({'sub_aspect': 1})

moll_fig_dimensions.update(generic_slice_fig_dimensions)
# these guys need to be joined not updated
molllonav_fig_dimensions =\
        dict(     list(generic_slice_fig_dimensions.items()) +\
                  list(molllonav_fig_dimensions.items())    )
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

    # shape to make geometric fields
    zero = np.zeros(np.array(np.shape(vals[..., 0])))

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
            cost = zero

        # first get root variable name and store any modifiers
        varname, deriv, primevar, sphvar = get_varprops(varname)

        # get sine/cotangent from cosine
        sint = np.sin(np.arccos(cost))
        cott = cost/sint
       
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
                the_slice += 2*eq.omega0
            elif varname[-1] == 'r': 
                the_slice += 2*eq.omega0*cost
            elif varname[-1] == 't':
                the_slice -= 2*eq.omega0*sint
            # for lambda and phi, there is no planetary vorticity 
            # contribution

        if 'pv' in varname: # need to divide by H
            # get rmin and rmax from the directory
            rmin, rmax = get_rminmax(dirname)
            the_slice /= compute_axial_H(rr, sint, rmin, rmax)

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
kw_plot_moll_or_ortho_default = dict({'clon': 0., 'clat': 20., 'shrinkage': 1., 'plotlonlines': True, 'lonvals': np.arange(0., 360., 30.), 'plotlatlines': True, 'latvals': np.arange(-60., 90., 30.), 'linewidth': 0.75*default_lw, 'plotboundary': True, 'ortho': False})
kw_plot_moll_or_ortho_default.update(kw_my_contourf_default)
kw_plot_moll_or_ortho_default['plotcontours'] = False

def plot_moll_or_ortho(field, costheta, fig, ax, **kw):
    kw = update_dict(kw_plot_moll_or_ortho_default, kw)
    find_bad_keys(kw_plot_moll_or_ortho_default, kw, 'plot_moll_or_otho')
    # change default plotcontours --> False in my_contourf
    kw_my_contourf = dotdict(kw_my_contourf_default.copy())
    kw_my_contourf.plotcontours = False
    if kw.ortho: # make colorbar a bit shorter
        kw_my_contourf.cbar_length_tol = 0.6

    print("my cf cbar_length_tol=", kw_my_contourf.cbar_length_tol)
    kw_my_contourf = update_dict(kw_my_contourf, kw)
    print("my cf cbar_length_tol=", kw_my_contourf.cbar_length_tol)
        
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
            ax.plot(xx_masked[imer, :], yy_masked[imer, :], 'k--', linewidth=linewidth)
    if kw.plotlatlines:
        for latval in kw.latvals:
            if latval == 0.: 
                linewidth = 2*kw.linewidth
            else:
                linewidth = kw.linewidth
            ilat = np.argmin(np.abs(lat - latval*np.pi/180.))
            ax.plot(xx_masked[:, ilat], yy_masked[:, ilat], 'k--', linewidth=linewidth)

    if kw.plotboundary:
        # Plot outer boundary
        psivals = np.linspace(0, 2*np.pi, 100)
        xvals, yvals = 2.*np.cos(psivals), np.sin(psivals)
        if kw.ortho:
            xvals /= 2.
        ax.plot(xvals, yvals, 'k', linewidth=1.5*kw.linewidth)

# equatorial slice plotting routine
kw_plot_eq_default = dict({'clon': 0., 'plotlonlines': True, 'lonvals': np.arange(0., 360., 60.), 'linewidth': 0.5*default_lw, 'plotboundary': True})
kw_plot_eq_default.update(kw_my_contourf_default)
kw_plot_eq_default['plotcontours'] = False

def plot_eq(field, rr, fig, ax, **kw):
    kw = update_dict(kw_plot_eq_default, kw)
    find_bad_keys(kw_plot_eq_default, kw, 'plot_eq')
    # change default plotcontours --> False in my_contourf
    tmp = kw_my_contourf_default.copy()
    tmp['plotcontours'] = False
    kw_my_contourf = update_dict(tmp, kw)
        
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
            ax.plot(xx[ilon, :], yy[ilon, :], 'k--', linewidth=linewidth)

    if kw.plotboundary:
        # plot boundaries
        ax.plot(xx[:, 0], yy[:, 0], 'k', linewidth=1.5*kw.linewidth)
        ax.plot(xx[:, -1], yy[:, -1], 'k', linewidth=1.5*kw.linewidth)

kw_plot_cutout_3d_default =\
        dotdict(dict({'r1': 'rmin', 'r2': 'rmax', 'dlon1': -30., 'dlon2': 60.,\
        'eq': True, 'varnames': 'vr', 'clon': 0., 'clat': 20., 't0': None, 'verbose': False, 'linewidth': default_lw, 'plotboundary': True,\
    'nocbar': True, 'rvals': [], 'lonvals': [], 'latvals': []}))
kw_plot_cutout_3d_default.update(kw_my_contourf_default)
kw_plot_cutout_3d_default.plotcontours = False

def plot_cutout_3d(dirname, fname, varname, fig, ax, **kw_in):
    kw = update_dict(kw_plot_cutout_3d_default, kw_in)
    kw_my_contourf = dotdict(kw_my_contourf_default.copy())
    kw_my_contourf.plotcontours = False
    kw_my_contourf = update_dict(kw_my_contourf, kw_in)
    find_bad_keys(kw_plot_cutout_3d_default, kw_in, 'plot_cutout_3d()')

    # get the Shell_Slice data
    ss = Shell_Slices(dirname + '/Shell_Slices/' + fname, '')
    mer = Meridional_Slices(dirname + '/Meridional_Slices/' + fname, '')
    if kw.eq:
        eq = Equatorial_Slices(dirname + '/Equatorial_Slices/' + fname, '')

    # unpack the actual data
    field_ss = get_slice(ss, varname, dirname=dirname)
    field_mer = get_slice(mer, varname, dirname=dirname)
    if kw.eq:
        field_eq = get_slice(eq, varname, dirname=dirname)

    # get gridinfo (do this individually for the particular instant,
    # in case resolution is changing)
    sliceinfo_ss = get_sliceinfo(dirname, 'Shell_Slices', fname=fname)
    sliceinfo_mer = get_sliceinfo(dirname, 'Meridional_Slices', fname=fname)
    gi = get_grid_info(dirname, ntheta=sliceinfo_mer.ntheta)

    # and unpack it
    ntheta = gi.nt
    nphi = gi.nphi
    lon1d = gi.phi - np.pi
    lat1d = gi.tt_lat*np.pi/180.

    # choose spherical cutout parameters

    # get radian versions of stuff
    clon, clat, dlon1, dlon2 = np.array([kw.clon, kw.clat, kw.dlon1, kw.dlon2])*np.pi/180.

    # these are desired longitudes
    lon1, lon2 = clon + dlon1, clon + dlon2

    # get these in range [-pi, pi)
    lon1, lon2, clon =\
            (np.array([lon1, lon2, clon]))%(2*np.pi) - np.pi

    # get indices and recompute values
    # clat and clon
    iclat = iclosest(lat1d, clat)
    clat = lat1d[iclat]
    iclon = iclosest(lon1d, clon)
    clon = lon1d[iclon]

    # shift the field_ss and field_eq so that the clon is in the 
    # just-right-of-center of the array (for [...) )
    # recall that ntheta = nphi/2, just-right-of-center is index ntheta
    field_ss = np.roll(field_ss, gi.nt - iclon, axis=0)
    if kw.eq:
        field_eq = np.roll(field_eq, gi.nt - iclon, axis=0)

    # get indices and recompute values for r1, r2, dlon1, dlon2
    r1, r2 = interpret_rvals(dirname, [kw.r1, kw.r2])
    rvals_ss = sliceinfo_ss.rvals
    ir1_ss, ir2_ss = inds_from_vals(rvals_ss, [r1, r2])
    r1, r2 = rvals_ss[[ir1_ss, ir2_ss]]
    ir1, ir2 = inds_from_vals(gi.rr, np.array([r1, r2]))
    beta = r1/r2 # this is the aspect ratio of the spherical cutout

    # need to shift available lons from Meridional Slices
    # to lie in range [-pi, pi)
    # apparently there is a mismatch between mer slices phi vals
    # and shell slice phivals - the indexing is different
    lonvals_mer = (sliceinfo_mer.lonvals + 0*np.pi)%(2*np.pi) - np.pi

    # now get the actual longitudes of the meridional slices we will use
    ilon1_mer, ilon2_mer = inds_from_vals(lonvals_mer, [lon1, lon2])
    # these are the actual longitudes we get
    lon1, lon2 = lonvals_mer[[ilon1_mer, ilon2_mer]]

    # now recompute dlon1, dlon2
    # in shifting everything around, they might not be in the correct range
    dlon1, dlon2 = lon1 - clon, lon2 - clon
    # make sure dlon1, dlon2 are negative and positive
    if dlon1 > 0.:
        dlon1 -= 2*np.pi
    if dlon2 < 0.:
        dlon2 += 2*np.pi

    # these are the near-final dlon1, dlon2 we end up 
    # now must check these are in the right range
    # if not, everything was for naught and we return to defaults!
    if dlon1 < -np.pi/2 or dlon1 > 0.:
        print("dlon1 = %.1f" %(dlon1*180./np.pi) + " is not allowed.")
        print("resetting to default dlon1 = -30.0")
        dlon1 = -np.pi/6
        ilon1 = iclosest(lonvals_mer, clon+dlon1)
        lon1 = lonvals_mer[ilon1]
        dlon1 = lon1 - clon
    if dlon2 < 0. or dlon2 > np.pi/2:
        print("dlon2 = %.1f" %(dlon2*180./np.pi) + " is not allowed.")
        print("resetting to default dlon2 = 60.0")
        dlon2 = np.pi/3
        ilon2 = iclosest(lonvals_mer, clon+dlon2)
        lon2 = lonvals_mer[ilon2]
        dlon2 = lon2 - clon

    # make sure dlon1, dlon2 are negative and positive (again)
    if dlon1 > 0.:
        dlon1 -= 2*np.pi
    if dlon2 < 0.:
        dlon2 += 2*np.pi

    # OK, finally done with getting correct longitudes

    # maybe print these in degrees
    if kw.verbose:
        # get degree versions of things
        clon_deg, clat_deg, lon1_deg, lon2_deg = np.array([clon, clat, lon1, lon2])*180./np.pi
        # print the parameters of the cutouts
        print(buff_line)
        print("PLOTTING SPHERICAL CUTOUTS!")
        print("clon =", clon_deg)
        print("clat =", clat_deg)
        print("r1 =", r1)
        print("r2 =", r2)
        print("(beta = " + str(beta) + ")")
        print("dlon1 =", lon1_deg - clon_deg)
        print("dlon2 =", lon2_deg - clon_deg)
        print("(lon1 =", lon1_deg, ")")
        print("(lon2 =", lon2_deg, ")")


    # Plot ORTHO 2 (the outer one)
    # get a "meshgrid" from 1D arrays
    lon2d, lat2d = np.meshgrid(lon1d, lat1d, indexing='ij')
    # get the field
    field = np.copy(field_ss[..., ir2_ss])
    if kw.nocbar:
        field /= np.std(field)
    # do ortho projection
    xx = np.cos(lat2d)*np.sin(lon2d)
    yy = np.cos(clat)*np.sin(lat2d) - np.sin(clat)*np.cos(lat2d)*np.cos(lon2d)
    # mask xx and yy to not plot INTERIOR triangle
    #don't plot the far side of the sphere
    cond1 = np.sin(clat)*np.sin(lat2d)+np.cos(clat)*np.cos(lat2d)*np.cos(lon2d) < 0
    # don't plot inside the meridional planes
    cond2 = (lon2d > dlon1) & (lon2d < dlon2)
    if kw.eq:
        # only don't plot inside the planes above the equator
        cond2 = cond2 & (lat2d > 0)
    field[cond1] = np.nan
    field[cond2] = np.nan

    # make the plot
    if kw.nocbar: # resolve this inconsistency
        kw_my_contourf.plotcbar = False
    my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)

    if kw.plotboundary:
        # plot the boundary
        nsvals = 1000 # s is "a parameter"

        # outer sphere
        svals = np.linspace(0, 2*np.pi, nsvals)
        xx, yy = np.cos(svals), np.sin(svals)
        if not kw.eq:
            # ignore the part of the boundary between the two meridional planes
            oldlon = np.arctan2(np.cos(svals), -np.sin(clat)*np.sin(svals))
            cond = (oldlon > dlon1) & (oldlon < dlon2)
            xx[cond] = np.nan
            yy[cond] = np.nan
        ax.plot(xx, yy, 'k-', linewidth=kw.linewidth)

    # PLOT MERIDIAN 1 (left one)
    # get the field
    field = np.copy(field_mer[ilon1_mer]).T
    if kw.nocbar: # saturate each radius separately
        for ir in range(gi.nr):
            field[ir, :] /= np.std(field[ir, :])
    # make mesh grid
    rad1d = np.copy(gi.rr)/r2
    rad2d, lat2d = np.meshgrid(rad1d, lat1d, indexing='ij')
    # make projection
    xx = rad2d * np.cos(lat2d) * np.sin(dlon1)
    yy = rad2d * (np.cos(clat)*np.sin(lat2d) - np.cos(dlon1)*np.sin(clat)*np.cos(lat2d))
    # mask outside the projection
    # don't plot outside the shell slices
    cond = (rad2d < beta)  | (rad2d > 1.)
    if kw.eq:
        # also don't plot below the equatorial plane
        cond = cond | (lat2d < 0)

    field[cond] = np.nan
    # make plot
    my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)

    # PLOT MERIDIAN 2 (right one)
    # get the field
    field = np.copy(field_mer[ilon2_mer, :, ::-1]).T
    if kw.nocbar: # saturate each radius separately
        for ir in range(gi.nr):
            field[ir, :] /= np.std(field[ir, :])
    # make mesh grid
    rad1d = np.copy(gi.rr[::-1])/r2
    rad2d, lat2d = np.meshgrid(rad1d, lat1d, indexing='ij')
    # make projection
    xx = rad2d * np.cos(lat2d) * np.sin(dlon2)
    yy = rad2d * (np.cos(clat)*np.sin(lat2d) - np.cos(dlon2)*np.sin(clat)*np.cos(lat2d))
    # mask outside the projection
    # don't plot outside the shell slices
    cond = (rad2d < beta)  | (rad2d > 1.)
    if kw.eq:
        # also don't plot below the equatorial plane
        cond = cond | (lat2d < 0)
    field[cond] = np.nan
    # make plot
    my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)

    if kw.plotboundary:
        # plot meridional boundaries
        # meridians 1 then 2, inner inner to outer
        rvals_extra = (make_array(kw.rvals)/r2).tolist()
        n_extra = len(rvals_extra)
        rvals = [beta, 1] + rvals_extra
        linestyles = 2*['-'] + n_extra*['--']
        nrvals = len(rvals)

        for dlonval in [dlon1, dlon2]:
            for irval in range(nrvals):
                if kw.eq:
                    svals = np.linspace(0., np.pi/2., nsvals) # is is latitude now
                else:
                    svals = np.linspace(-np.pi/2, np.pi/2., nsvals)
                rval = rvals[irval]

                # get the xx and yy points
                xx = rval*np.cos(svals)*np.sin(dlonval)
                yy = rval*(np.cos(clat)*np.sin(svals) - np.sin(clat)*np.cos(svals)*np.cos(dlonval))
            
                # don't plot anything behind the inner sphere (i.e., same condition for ortho masks)
                if rval < beta/np.cos(clat):
                    cond = np.sin(clat)*np.sin(svals)+np.cos(clat)*np.cos(svals)*np.cos(dlonval) < 0
                    xx[cond] = np.nan
                    yy[cond] = np.nan

                ax.plot(xx, yy, 'k', linewidth=kw.linewidth, linestyle=linestyles[irval])

    # Plot the ORTHO 1 (the inner one)
    # get a "meshgrid" from 1D arrays
    lon2d, lat2d = np.meshgrid(lon1d, lat1d, indexing='ij')
    # get the field
    field = np.copy(field_ss[..., ir1_ss])
    if kw.nocbar: # "saturate separately"
        field /= np.std(field)

    # do ortho projection
    xx = beta*np.cos(lat2d)*np.sin(lon2d)
    yy = beta*(np.cos(clat)*np.sin(lat2d) - np.sin(clat)*np.cos(lat2d)*np.cos(lon2d))
    # mask xx and yy such that only the INNER TRIANGLE is visible
    #don't plot the far side of the sphere
    cond1 = np.sin(clat)*np.sin(lat2d)+np.cos(clat)*np.cos(lat2d)*np.cos(lon2d) < 0
     # don't plot outside the meridional planes
    cond2 = (lon2d < dlon1) | (lon2d > dlon2)
    if kw.eq:
        # also don't plot below the equatorial plane
        cond2 = cond2 | (lat2d < 0)
    field[cond1] = np.nan
    field[cond2] = np.nan
    if kw.eq:
        # only don't plot inside the planes above the equatorial plane
        cond2 = cond2 & (lat2d > 0)
    my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)

    if not kw.eq and kw.plotboundary:
        # plot the boundary
        # inner sphere
        svals = np.linspace(0, 2*np.pi, nsvals)
        xx, yy = beta*np.cos(svals), beta*np.sin(svals)
        if not kw.eq:
            # ignore the part of the boundary outside the two meridional planes
            oldlon = np.arctan2(np.cos(svals), -np.sin(clat)*np.sin(svals))
            cond = (oldlon < dlon1) | (oldlon > dlon2)
            xx[cond] = np.nan
            yy[cond] = np.nan
        ax.plot(xx, yy, 'k-', linewidth=kw.linewidth)

    # PLOT EQUATORIAL SLICE
    if kw.eq:
        # get field
        field = np.copy(field_eq)
        if kw.nocbar: # saturate each radius separately
            for ir in range(gi.nr):
                field[:, ir] /= np.std(field[:, ir])
        # make mesh grid
        rad1d = np.copy(gi.rr)/r2
        lon2d, rad2d = np.meshgrid(lon1d, rad1d, indexing='ij')
        # make projection
        xx = rad2d*np.sin(lon2d)
        yy = -rad2d*np.sin(clat)*np.cos(lon2d)
        # mask outside the projection
        cond = (lon2d < dlon1) | (lon2d > dlon2) | (rad2d < beta)  | (rad2d > 1.)
        field[cond] = np.nan
        # make plot
        my_contourf(xx, yy, field, fig, ax, **kw_my_contourf)

    # don't obscure the boundary
    lilbit = 0.01
    ax.set_xlim(-1-lilbit, 1 + lilbit)
    ax.set_ylim(-1-lilbit, 1+lilbit)

    # then finish outer meridians with dashed lines
    if kw.plotboundary:
        count = 0
        for dlonval in [dlon1, dlon2, dlon2-np.pi, dlon1+np.pi]:
            if count < 2:
                svals = np.linspace(-np.pi/2., 0., nsvals)
            else:
                svals = np.linspace(0., np.pi/2., nsvals)
            xx, yy = np.cos(svals)*np.sin(dlonval), np.cos(clat)*np.sin(svals) - np.sin(clat)*np.cos(svals)*np.cos(dlonval)
            cond = np.sin(clat)*np.sin(svals)+np.cos(clat)*np.cos(svals)*np.cos(dlonval) < 0
            xx[cond] = np.nan
            yy[cond] = np.nan
            ax.plot(xx, yy, 'k--', linewidth=kw.linewidth)
            count += 1

        # and the "far" meridians
        #svals = np.linspace(0., np.pi/2., nsvals)
        #for dlonval in [dlon2-np.pi, dlon1+np.pi]:
        #    xx, yy = np.cos(svals)*np.sin(dlonval), np.cos(clat)*np.sin(svals) - np.sin(clat)*np.cos(svals)*np.cos(dlonval)
        #    ax.plot(xx, yy, 'k--', linewidth=kw.linewidth)

        # plot equatorial plane boundaries

        # constant longitude lines in equatorial plane
        if kw.eq:
            linestyle = '-'
        else:
            linestyle = '--'
        svals = np.linspace(beta, 1, nsvals)
        for lonval in [dlon1, dlon2]:
            ax.plot(svals*np.sin(lonval), -svals*np.sin(clat)*np.cos(lonval),\
                    'k', linewidth=kw.linewidth, linestyle=linestyle)

        # equatorial slice, inner to outer
        if kw.eq:
            irvals = np.arange(nrvals)
            linestyles_loc = np.copy(linestyles)
        else: # just plot outer and inner boundaries
            irvals = [0]
            linestyles_loc = ['--']
        svals = np.linspace(dlon1, dlon2, nsvals)
        for irval in irvals:
            rval = rvals[irval]
            ax.plot(rval*np.sin(svals), -rval*np.sin(clat)*np.cos(svals), 'k',\
                    linewidth=kw.linewidth, linestyle=linestyles_loc[irval])
            
        # finish the equator with a dashed line 
        # (this happens no matter if the equatorial plane is there or not)
        svals = np.linspace(dlon1, dlon2, nsvals)
        xx, yy = np.sin(lon1d), -np.sin(clat)*np.cos(lon1d)
        # mask outside the projection
        cond1 = np.cos(clat)*np.cos(lon1d) < 0
        cond2 = (lon1d > dlon1) & (lon1d < dlon2)
        xx[cond1] = np.nan; xx[cond2] = np.nan
        yy[cond1] = np.nan; yy[cond2] = np.nan
        ax.plot(xx, yy, 'k--', linewidth=kw.linewidth)

        # upper rotation axis
        svals = np.linspace(beta, 1, nsvals)
        ax.plot(0*svals, svals*np.cos(clat), 'k-', linewidth=kw.linewidth)
        if not kw.eq: # lower rotation axis too
            svals = np.linspace(-1, -beta/np.cos(clat), nsvals)
            ax.plot(0*svals, svals*np.cos(clat), 'k-', linewidth=kw.linewidth)
