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

def plot_ortho(field_orig, radius, costheta, fig=None, ax=None, ir=0,\
        minmax=None, clon=0, clat=20, posdef=False, logscale=False,\
        lw_scaling=1., plot_cbar=True, cbar_fs=10,\
        symlog=False, linscale=None, linthresh=None, cmap=None,\
        bold_patch=None, thickcenter=True, cbar_scaling=1.):
    
    if logscale:
        posdef = True
        
    # Shouldn't have to do this but Python is stupid with arrays ...
    field = np.copy(field_orig)
    
    # Get geometric parameters
    ro = np.max(radius)
    rloc = radius[ir]
    
    # Get latitude and longitude grid from the costheta array
    theta = np.arccos(costheta)
    lats = 90. - theta*180./np.pi
    ntheta = len(theta)
    nphi = 2*ntheta
    dlon = 360.0/(nphi - 1.) # this eliminates a whitespace slice
    # between the first and last elements of the shifted data array
    # Technically this is a distortion, but only 1 part in nphi!
    lons = dlon*np.arange(nphi)  # In rayleigh, lons[0] = 0.
    # Make 2d (meshed) grid of longitude/latitude
    llon, llat = np.meshgrid(lons, lats, indexing='ij')
    
    # Essence of cartopy is transformations:
    # This one will be from PlateCarree (simple longitude/latitude equally 
    # spaced) to orthographic
    pc = crs.PlateCarree()
    ortho = crs.Orthographic(central_latitude=clat, central_longitude=180.)
    
	# "Roll" the original field so that so that the desired clon falls
	# on 180 degrees, which is the ~center of the lons array
    difflon = clon - 180.
    iphi_shift = -int(difflon/360.*nphi)
    field = np.roll(field, iphi_shift, axis=0)

    # Get default bounds if not specified
    if minmax is None:
        minmax = get_satvals(field, posdef=posdef, logscale=logscale,\
                        symlog=symlog)

    # Get the exponent to use for scientific notation
    if not (logscale or symlog):
        maxabs = max(np.abs(minmax[0]), np.abs(minmax[1]))
        maxabs_exp = int(np.floor(np.log10(maxabs)))
        divisor = 10**maxabs_exp
        
        # Normalize field by divisor
        field /= divisor
        minmax = minmax[0]/divisor, minmax[1]/divisor
        # Can't reassign tuples element-wise for some reason

    # Saturate the array (otherwise creates white space with contourf)
    saturate_array(field, minmax[0], minmax[1])     
         
    # Get the orthographic projection coordinates of llon/llat by 
    # converting between PlateCarree (lat/lon) and orthographic
    x_withnan = np.zeros((nphi, ntheta)) 
    y_withnan = np.zeros((nphi, ntheta)) 
    for iphi in range(nphi):
        for itheta in range(ntheta):
            lon_loc, lat_loc = llon[iphi, itheta], llat[iphi, itheta]
            x_withnan[iphi, itheta], y_withnan[iphi, itheta] =\
                ortho.transform_point(lon_loc, lat_loc, pc)

    # Must deal with "nan" coordinates since part of sphere is invisible
    x_unshrunk, y_unshrunk = deal_with_nans(x_withnan, y_withnan)
    coord_radius = np.sqrt(np.max(x_unshrunk**2. + y_unshrunk**2.))
    x_unshrunk /= coord_radius
    y_unshrunk /= coord_radius
    # Given cartopy's convention, coord_radius will generally be
    # ~6.37 x 10^6 (Earth radius in m), however can vary somewhat. It is
    # important to normalize the latitude and longitude lines by this 
    # amount (and the shrink factor) in order to display them properly!

    # Shrink the normalized coordinates to give depth perception
    shrink_factor = rloc/ro
    x = x_unshrunk*shrink_factor #+ shrink_distance
    y = y_unshrunk*shrink_factor #+ shrink_distance

    # Create a default fig/ax pair, if calling routine didn't specify
    # them
    figwasNone = False
    if fig is None or ax is None:
        fig, ax = default_axes_1by1()
        figwasNone = True

    ax.set_xlim((-1.01, 1.01)) # deal with annoying whitespace cutoff issue
    ax.set_ylim((-1.01, 1.01))
    ax.axis('off') # get rid of x/y axis coordinates

    # Set the colormap if unspecified
    if cmap is None:
        if logscale:
            cmap = 'Greys'
        elif posdef:
            cmap = 'plasma'
        elif symlog:
            cmap = 'RdYlBu_r'
        else:
            cmap = 'RdYlBu_r'

    # Plot the orthographic projection
    if logscale:
        log_min, log_max = np.log10(minmax[0]), np.log10(minmax[1])
        levels = np.logspace(log_min, log_max, 150)
        im = ax.contourf(x, y, field, cmap=cmap,\
            norm=colors.LogNorm(vmin=minmax[0], vmax=minmax[1]),\
            levels=levels)  
    elif posdef:
        levels = np.linspace(minmax[0], minmax[1])
        im = ax.contourf(x, y, field, cmap=cmap, levels=levels)
    elif symlog:
        linthresh_default, linscale_default =\
            get_symlog_params(field, field_max=minmax[1])
        if linthresh is None:
            linthresh = linthresh_default
        if linscale is None:
            linscale = linscale_default
        log_thresh = np.log10(linthresh)
        log_max = np.log10(minmax[1])
        nlevs_per_interval = 100
        levels_neg = -np.logspace(log_max, log_thresh,\
                nlevs_per_interval,\
                endpoint=False)
        levels_mid = np.linspace(-linthresh, linthresh,\
                nlevs_per_interval, endpoint=False)
        levels_pos = np.logspace(log_thresh, log_max,\
                nlevs_per_interval)
        levels = np.hstack((levels_neg, levels_mid, levels_pos))
        im = ax.contourf(x, y, field, cmap=cmap,\
            norm=colors.SymLogNorm(linthresh=linthresh,\
            linscale=linscale, vmin=minmax[0], vmax=minmax[1]),\
            levels=levels)
    else:
        im = ax.contourf(x, y, field, cmap=cmap,\
                levels=np.linspace(minmax[0], minmax[1], 150))        
       
    # Draw parallels and meridians, evenly spaced by 30 degrees
    default_lw = 0.5*lw_scaling # default linewidth bit thinner
    parallels = np.arange(-60., 90., 30.)
    
    # Make sure the plotted meridians take into account the shift
    # in the phi-coordinate of the data: data that was at phi = 0. --> 
    # data at phi = -difflon,
    # i.e, our phi = 0. line should appear at -difflon
    if lonvals is None:
        lonvals = np.arange(-difflon, -difflon + 360., 30.)

    npoints = 100
    for lonval in lonvals:
        if meridian == -difflon and thickcenter: 
            lw = 1.3*lw_scaling # make central longitude thicker

        else:
            lw = default_lw
            
        lonvals = meridian + np.zeros(npoints)
        latvals = np.linspace(-90, 90, npoints)
        linex_withnans, liney_withnans = np.zeros(npoints),\
                np.zeros(npoints)
        for i in range(npoints):
            linex_withnans[i], liney_withnans[i] =\
                ortho.transform_point(lonvals[i], latvals[i], pc)
        linex, liney =\
            deal_with_nans(linex_withnans.reshape((npoints, 1)),\
            liney_withnans.reshape((npoints, 1)))
        # normalize by coord_radius ~ Earth radius:
        linex = linex[:, 0]/coord_radius
        liney = liney[:, 0]/coord_radius

        # Apply the shrink factor to give depth perception
        linex *= shrink_factor
        liney *= shrink_factor

        # Plot a thin black line
        ax.plot(linex, liney, 'k', linewidth=lw)

    for parallel in parallels:
        if parallel == 0. and thickcenter:
                lw = 1.3*lw_scaling # make equator thicker
        else:
            lw = default_lw        
        lonvals = np.linspace(0., 360., npoints)
        latvals = np.zeros(npoints) + parallel
        linex_withnans, liney_withnans =\
            np.zeros(npoints), np.zeros(npoints)
        for i in range(npoints):
            linex_withnans[i], liney_withnans[i] =\
                ortho.transform_point(lonvals[i], latvals[i], pc)
        linex, liney =\
            deal_with_nans(linex_withnans.reshape((npoints, 1)),\
            liney_withnans.reshape((npoints, 1)))
        # normalize by coord_radius ~ Earth radius:
        linex = linex[:, 0]/coord_radius
        liney = liney[:, 0]/coord_radius

        # Apply the shrink factor to give depth perception
        linex *= shrink_factor
        liney *= shrink_factor

        # Plot a thin black line
        ax.plot(linex, liney, 'k', linewidth=lw)   

    # Possibly highlight a patch of the lat/lon grid
    if not bold_patch is None:
        lon0, lat0 = bold_patch # lower-left corner of that patch as l
        lonvals = np.linspace(lon0, lon0 + 30., npoints)
        latvals = np.linspace(lat0, lat0 + 30., npoints)
        for parallel in [lat0, lat0 + 30.]:
            linex_withnans, liney_withnans =\
                np.zeros(npoints), np.zeros(npoints)
            for i in range(npoints):
                linex_withnans[i], liney_withnans[i] =\
                    ortho.transform_point(lonvals[i], parallel, pc)
            linex, liney =\
                deal_with_nans(linex_withnans.reshape((npoints, 1)),\
                liney_withnans.reshape((npoints, 1)))
            # normalize by coord_radius ~ Earth radius:
            linex = linex[:, 0]/coord_radius
            liney = liney[:, 0]/coord_radius

            # Apply the shrink factor to give depth perception
            linex *= shrink_factor
            liney *= shrink_factor

            # Plot a thick black line
            lw = 1.5*lw_scaling # make central longitude thicker
            ax.plot(linex, liney, 'k', linewidth=lw)   

        for lon in [lon0, lon0 + 30.]:
            linex_withnans, liney_withnans =\
                np.zeros(npoints), np.zeros(npoints)
            for i in range(npoints):
                linex_withnans[i], liney_withnans[i] =\
                    ortho.transform_point(lon, latvals[i], pc)
            linex, liney =\
                deal_with_nans(linex_withnans.reshape((npoints, 1)),\
                liney_withnans.reshape((npoints, 1)))
            # normalize by coord_radius ~ Earth radius:
            linex = linex[:, 0]/coord_radius
            liney = liney[:, 0]/coord_radius

            # Apply the shrink factor to give depth perception
            linex *= shrink_factor
            liney *= shrink_factor

            # Plot a thick black line
            lw = 1.5*lw_scaling # make central longitude thicker
            ax.plot(linex, liney, 'k', linewidth=lw)   

    # Set up color bar
    if plot_cbar:
        ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
        ax_delta_x = ax_xmax - ax_xmin
        ax_delta_y = ax_ymax - ax_ymin  
    
        cbar_aspect = 1./20.
        fig_aspect = ax_delta_x/ax_delta_y # assumes subplot aspect ratio is 1
            # otherwise, must multiply by proper subplot aspect ratio (=0.5 for
            # Mollweide))
        cbar_width = 0.5*ax_delta_x*cbar_scaling # make cbar half as long as plot is wide
        cbar_height = cbar_width*cbar_aspect/fig_aspect
        cbar_bottom = ax_ymin - 2.5*cbar_height
        cbar_left = ax_xmin + 0.5*ax_delta_x - 0.5*cbar_width
    
        cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))
        cbar = plt.colorbar(im, cax=cax, orientation='horizontal')    
        
        if logscale:
            locator = ticker.LogLocator(subs='all')
            cbar.set_ticks(locator)
            cbar_units = ' ' + units #texunits.get(varname, 'cgs')
        elif posdef:
            cbar_units = ' ' + (r'$\times10^{%i}$' %maxabs_exp) + ' ' + units#\
                #texunits.get(varname, 'cgs')
            cbar.set_ticks([minmax[0], minmax[1]])
            cbar.set_ticklabels(['%1.2f' %minmax[0], '%1.2f' %minmax[1]])
        elif symlog:
            cbar_units = ' ' + units #texunits.get(varname, 'cgs')
            nlin = 5
            nlog = 6
            lin_ticks = np.linspace(-linthresh, linthresh, nlin)
            log_ticks1 = np.linspace(minmax[0], -linthresh, nlog,\
                    endpoint=False)
            log_ticks2 = -log_ticks1[::-1]
            ticks = np.hstack((log_ticks1, lin_ticks, log_ticks2))
            nticks = nlin + 2*nlog
            cbar.set_ticks(ticks)
            ticklabels = []
            for i in range(nticks):
                ticklabels.append(r'')
            ticklabels[0] = sci_format(minmax[0])
            ticklabels[nlog] = sci_format(-linthresh)
            ticklabels[nticks//2] = r'$0$'
            ticklabels[nlog + nlin - 1] = sci_format(linthresh)
            ticklabels[nticks - 1] = sci_format(minmax[1])
            cbar.set_ticklabels(ticklabels)
#            cbar.set_ticks([-minmax[1], -linthresh, 0, linthresh,\
#                    minmax[1]])
#            cbar.set_ticklabels([sci_format(-minmax[1]),\
#                    sci_format(-linthresh), '0', sci_format(linthresh),\
#                    sci_format(minmax[1])])
    #            cax.minorticks_on()
        else:
            cbar_units = ' ' + (r'$\times10^{%i}$' %maxabs_exp) + units#\
                    #' ' + texunits.get(varname, 'cgs')
            cbar.set_ticks([minmax[0], 0, minmax[1]])
            cbar.set_ticklabels(['%1.2f' %minmax[0], '0', '%1.2f'\
                    %minmax[1]])
        # Title the colorbar based on the field's units
        fig.text(cbar_left + cbar_width, cbar_bottom + 0.5*cbar_height,\
                 cbar_units, verticalalignment='center',\
                 fontsize=cbar_fs)     
        plt.sca(cax)
        plt.xticks(fontsize=cbar_fs) # change the fontsize 
                            # on the tick labels too

    # Plot outer boundary
    psivals = np.linspace(0, 2*np.pi, 500)
    xvals, yvals = np.cos(psivals), np.sin(psivals)
    ax.plot(xvals, yvals, 'k', linewidth=1.3*lw_scaling)

    if figwasNone: # user probably called plot_ortho from the python 
        # command line, wanting to view the projection immediately
        plt.show()
    del field # free up memory
    return im
