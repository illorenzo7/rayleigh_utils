import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib import ticker, colors
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from cartopy import crs
from varprops import texunits
from common import get_satvals, saturate_array, sci_format, get_symlog_params
from get_parameter import get_parameter

def deal_with_nans(x, y):
    """ Nandles the NaN in x and y produced by an orthographic projection
    by making quadrilateral corners on the "other side" of the projection 
    (which the NaNs correspond to) overlap with adjacent quadrilateral
    corners--so that pcolormesh and contourf plot no data from the "far 
    side" of the sphere. 
    Note that for an orthographic projection,
    x (East-West map coordinate) ~ longitude and y (North-South map
     coordinate) ~ latitude.
    Also in plot_ortho, the 0th axis of the arrays corresponds to longitude,
    and the 1st axis corresponds to latitude.
    =========================================
    NOTE: deal_with_nans MUST be passed arrays in which the nans occur in
    the first and last parts of each column (range of x coords, or 
    longitudes)--i.e., roll the array so that your desired clon appears 
    in the middle of the array BEFORE transforming coordinates with 
    cartopy. If the nans occupy the middle of the column, pcolormesh
    and contourf will do bad things.
    """
    x_dealt, y_dealt = np.copy(x), np.copy(y)

    found_min_iy = False
    found_max_iy = False
    
    nx, ny = np.shape(x)
    min_iy, max_iy = 0, ny 
    # assume by default no max and min iy 
    # ("latitude") values    
    for iy in range(ny):
        xvals = np.copy(x_dealt[:, iy])
        good = np.where(np.isfinite(xvals))[0]
        if good.size == 0:
            if found_min_iy and not found_max_iy:
                max_iy = iy - 1
                found_max_iy = True
            continue
        elif good.size == nx: # no nans, so nothing to do, so move on
            # if this is the first column with "good" (non-nan) elements,
            # you've found the min. iy!
            if not found_min_iy:
                min_iy = iy
                found_min_iy = True
            continue
        else: # There are some nans to deal with
            # These had better occur in blocks at the beginning and end
            # of column iy. You also may have found min_iy!
            if not found_min_iy:
                min_iy = iy
                found_min_iy = True            
            
            # Find the nans and replace them with their closest x-coord
            # (~longitude)
            ix_min, ix_max = good[0], good[-1]
            x_dealt[:ix_min, iy] = x_dealt[ix_min, iy]
            x_dealt[ix_max + 1:, iy] = x_dealt[ix_max, iy]
            y_dealt[:ix_min, iy] = y_dealt[ix_min, iy]
            y_dealt[ix_max + 1:, iy] = y_dealt[ix_max, iy]
    
    for ix in range(nx):
        # Replace the nans before min_iy
        x_dealt[ix, :min_iy] = x_dealt[ix, min_iy]
        y_dealt[ix, :min_iy] = y_dealt[ix, min_iy]
        # ...and after max_iy, if applicable
        if max_iy < ny - 1:
            x_dealt[ix, max_iy + 1:] = x_dealt[ix, max_iy]
            y_dealt[ix, max_iy + 1:] = y_dealt[ix, max_iy]            
       
    return x_dealt, y_dealt


def rbounds(dirname):
    try:
        rmin = get_parameter(dirname, 'rmin')
        rmax = get_parameter(dirname, 'rmax')
#    if rmin == 100 or rmax == 100: # get_parameter must have failed
    except:
        dom_bounds = get_parameter(dirname, 'domain_bounds')
        rmin, rmax = dom_bounds[0], dom_bounds[-1]
    return rmin, rmax

def axis_range(ax): # gets subplot coordinates on a figure in "normalized"
        # coordinates
    pos = plt.get(ax, 'position')
    bottom_left = pos.p0
    top_right = pos.p1
    xmin, xmax = bottom_left[0], top_right[0]
    ymin, ymax = bottom_left[1], top_right[1]
    
    return xmin, xmax, ymin, ymax
    
def default_axes_1by1():
    # Create plot
    subplot_width_inches = 7.5
    subplot_height_inches = 7.5
    margin_inches = 1./8.
    margin_bottom_inches = 3./4.

    fig_width_inches = subplot_width_inches + 2.*margin_inches
    fig_height_inches = subplot_height_inches + margin_inches +\
            margin_bottom_inches

    margin_x = margin_inches/fig_width_inches
    margin_bottom = margin_bottom_inches/fig_height_inches
    subplot_width = subplot_width_inches/fig_width_inches
    subplot_height = subplot_height_inches/fig_height_inches

    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax = fig.add_axes((margin_x, margin_bottom, subplot_width,\
            subplot_height))
    return fig, ax

def default_axes_2by1():
    # Create plot
    subplot_width_inches = 10.
    subplot_height_inches = 5.
    margin_inches = 1./8.
    margin_bottom_inches = 3./4. # leave room for colorbar

    fig_width_inches = subplot_width_inches + 2.*margin_inches
    fig_height_inches = subplot_height_inches + margin_inches +\
            margin_bottom_inches

    margin_x = margin_inches/fig_width_inches
    margin_bottom = margin_bottom_inches/fig_height_inches
    subplot_width = subplot_width_inches/fig_width_inches
    subplot_height = subplot_height_inches/fig_height_inches

    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax = fig.add_axes((margin_x, margin_bottom, subplot_width,\
            subplot_height))
    return fig, ax

def plot_ortho(field_orig, radius, costheta, fig=None, ax=None, ir=0,\
        minmax=None, clon=0, clat=20, posdef=False, logscale=False,\
        varname='vr', lw_scaling=1., plot_cbar=True, cbar_fs=10,\
        symlog=False, linscale=None, linthresh=None):
    
    if 'sq' in varname or logscale:
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
        
    # Plot the orthographic projection
    if logscale:
        log_min, log_max = np.log10(minmax[0]), np.log10(minmax[1])
        levels = np.logspace(log_min, log_max, 150)
        im = ax.contourf(x, y, field, cmap='Greys',\
            norm=colors.LogNorm(vmin=minmax[0], vmax=minmax[1]),\
            levels=levels)  
    elif posdef:
        levels = np.linspace(minmax[0], minmax[1])
        im = ax.contourf(x, y, field, cmap='plasma', levels=levels)
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
        im = ax.contourf(x, y, field, cmap='RdYlBu_r',\
            norm=colors.SymLogNorm(linthresh=linthresh,\
            linscale=linscale, vmin=minmax[0], vmax=minmax[1]),\
            levels=levels)
    else:
        im = ax.contourf(x, y, field, cmap='RdYlBu_r',\
                levels=np.linspace(minmax[0], minmax[1], 150))        
       
    # Draw parallels and meridians, evenly spaced by 30 degrees
    default_lw = 0.5*lw_scaling # default linewidth bit thinner
    parallels = np.arange(-60., 90., 30.)
    
    # Make sure the plotted meridians take into account the shift
    # in the phi-coordinate of the data: data that was at phi = 0. --> 
    # data at phi = -difflon,
    # i.e, our phi = 0. line should appear at -difflon
    meridians = np.arange(-difflon, -difflon + 360., 30.)

    npoints = 100
    for meridian in meridians:
        if meridian == -difflon: 
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
        if parallel == 0.: 
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

    # Set up color bar
    if plot_cbar:
        ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
        ax_delta_x = ax_xmax - ax_xmin
        ax_delta_y = ax_ymax - ax_ymin  
    
        cbar_aspect = 1./20.
        fig_aspect = ax_delta_x/ax_delta_y # assumes subplot aspect ratio is 1
            # otherwise, must multiply by proper subplot aspect ratio (=0.5 for
            # Mollweide))
        cbar_width = 0.5*ax_delta_x # make cbar half as long as plot is wide
        cbar_height = cbar_width*cbar_aspect/fig_aspect
        cbar_bottom = ax_ymin - 2.5*cbar_height
        cbar_left = ax_xmin + 0.5*ax_delta_x - 0.5*cbar_width
    
        cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))
        cbar = plt.colorbar(im, cax=cax, orientation='horizontal')    
        
        if logscale:
            locator = ticker.LogLocator(subs='all')
            cbar.set_ticks(locator)
            cbar_units = ' ' + texunits[varname]
        elif posdef:
            cbar_units = ' ' + (r'$\times10^{%i}$' %maxabs_exp) + ' ' +\
                texunits[varname]
            cbar.set_ticks([minmax[0], minmax[1]])
            cbar.set_ticklabels(['%1.1f' %minmax[0], '%1.1f' %minmax[1]])
        elif symlog:
            cbar_units = ' ' + texunits[varname]
            cbar.set_ticks([-minmax[1], -linthresh, 0, linthresh,\
                    minmax[1]])
            cbar.set_ticklabels([sci_format(-minmax[1]),\
                    sci_format(-linthresh), '0', sci_format(linthresh),\
                    sci_format(minmax[1])])
    #            cax.minorticks_on()
        else:
            cbar_units = ' ' + (r'$\times10^{%i}$' %maxabs_exp) +\
                    ' ' + texunits[varname]
            cbar.set_ticks([minmax[0], 0, minmax[1]])
            cbar.set_ticklabels(['%1.1f' %minmax[0], '0', '%1.1f'\
                    %minmax[1]])
        # Title the colorbar based on the field's units
        fig.text(cbar_left + cbar_width, cbar_bottom + 0.5*cbar_height,\
                 cbar_units, verticalalignment='center', **csfont,\
                 fontsize=cbar_fs)     

    # Plot outer boundary
    psivals = np.linspace(0, 2*np.pi, 500)
    xvals, yvals = np.cos(psivals), np.sin(psivals)
    ax.plot(xvals, yvals, 'k', linewidth=1.3*lw_scaling)

    if figwasNone: # user probably called plot_ortho from the python 
        # command line, wanting to view the projection immediately
        plt.show()

def plot_moll(field_orig, costheta, fig=None, ax=None, minmax=None,\
        clon=0., posdef=False, logscale=False, symlog=False, varname='vr',\
        lw_scaling=1., plot_cbar=True, cbar_fs=10, linscale=None,\
        linthresh=None): 
    
    if 'sq' in varname or logscale:
        posdef = True
        
    # Shouldn't have to do this but Python is stupid with arrays ...
    field = np.copy(field_orig)    

    # Set tick label sizes (for colorbar)
    mpl.rcParams['xtick.labelsize'] = cbar_fs
    mpl.rcParams['ytick.labelsize'] = cbar_fs
    
    # Get latitude and longitude grid from the costheta array
    theta = np.arccos(costheta)
    lats = 90. - theta*180./np.pi
    ntheta = len(theta)
    nphi = 2*ntheta
    dlon = 360.0/nphi     
    lons = dlon*np.arange(nphi)  # In rayleigh, lons[0] = 0.
    # Make 2d (meshed) grid of longitude/latitude
    llon, llat = np.meshgrid(lons, lats, indexing='ij')    
    
    # Essence of cartopy is transformations:
    # This one will be from PlateCarree (simple longitude/latitude equally 
    # spaced) to Mollweide
    pc = crs.PlateCarree()
    moll = crs.Mollweide(central_longitude=180.)

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
    
    # Saturate the array (otherwise contourf will show white areas)
    saturate_array(field, minmax[0], minmax[1])    
            
    # Get the Mollweide projection coordinates of llon/llat by converting
    # between PlateCarree (lat/lon) and Mollweide
    points_moll = moll.transform_points(pc, llon, llat) 
        # must see why this doesn't work for othographic proj.
    x = points_moll[:, :, 0]
    y = points_moll[:, :, 1]

    # "normalize coordinates":
    radius = max(np.max(np.abs(x))/2., np.max(np.abs(y)))
    x /= radius # now falls in [-2, 2]
    y /= radius # now falls in [-1, 1]

    # Create a default fig/ax pair, if calling routine didn't specify
    # them
    figwasNone = False
    if fig is None and ax is None:
        fig, ax = default_axes_2by1()
        figwasNone = True

    ax.set_xlim((-2.02, 2.02)) # deal with annoying whitespace cutoff issue
    ax.set_ylim((-1.01, 1.01))
    ax.axis('off') # get rid of x/y axis coordinates
            
    # Make the Mollweide projection
    if logscale:
        log_min, log_max = np.log10(minmax[0]), np.log10(minmax[1])
        levels = np.logspace(log_min, log_max, 150)
        im = ax.contourf(x, y, field, cmap='Greys',\
            norm=colors.LogNorm(vmin=minmax[0], vmax=minmax[1]),\
            levels=levels)  
    elif posdef:
        levels = np.linspace(minmax[0], minmax[1])
        im = ax.contourf(x, y, field, cmap='plasma', levels=levels)
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
        im = ax.contourf(x, y, field, cmap='RdYlBu_r',\
            norm=colors.SymLogNorm(linthresh=linthresh,\
            linscale=linscale, vmin=minmax[0], vmax=minmax[1]),\
            levels=levels)
    else:
        im = ax.contourf(x, y, field, cmap='RdYlBu_r',\
                levels=np.linspace(minmax[0], minmax[1], 150))
      
    # Draw parallels and meridians, evenly spaced by 30 degrees
    default_lw = 0.5*lw_scaling # default linewidth bit thinner
    parallels = np.arange(-60., 90., 30.)
    
    # Make sure the plotted meridians take into account the shift
    # in the phi-coordinate of the data: data that was at phi = 0. --> 
    # data at phi = -difflon,
    # i.e, our phi = 0. line should appear at -difflon
    meridians = np.arange(-difflon, -difflon + 360., 30.)
    
    npoints = 100
    for meridian in meridians:
        if meridian == -difflon: 
            lw = 1.3*lw_scaling # make central longitude thicker
        else:
            lw = default_lw
            
        lonvals = meridian + np.zeros(npoints)
        latvals = np.linspace(-90, 90, npoints)
        linepoints = moll.transform_points(pc, lonvals, latvals)
        linex, liney = linepoints[:, 0], linepoints[:, 1]
        linex /= radius
        liney /= radius
        
        ax.plot(linex, liney, 'k', linewidth=lw)
    
    for parallel in parallels:
        if parallel == 0.: 
            lw = 1.3*lw_scaling # make equator thicker
        else:
            lw = default_lw        
        lonvals = np.linspace(0., 359.99, npoints)
        # lon = 0 = 360 is treated as a negative-x point for clon=180
        # so just below 360 is the positive-most x point
        latvals = np.zeros(npoints) + parallel
        linepoints = moll.transform_points(pc, lonvals, latvals)
        linex, liney = linepoints[:, 0], linepoints[:, 1]
        linex /= radius
        liney /= radius
        
        ax.plot(linex, liney, 'k', linewidth=lw)

    # Set up color bar
    if plot_cbar:
        ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
        ax_delta_x = ax_xmax - ax_xmin
        if symlog:
            cbar_aspect = 1./40.
            cbar_width = 0.5*ax_delta_x # make cbar one half 
                            # as long as plot is wide            
        else:
            cbar_aspect = 1./20.
            cbar_width = 0.25*ax_delta_x # make cbar one quarter
                            # as long as plot is wide            
        fig_width_inches, fig_height_inches = fig.get_size_inches()
        fig_aspect = fig_height_inches/fig_width_inches
        cbar_height = cbar_width*cbar_aspect/fig_aspect
        cbar_bottom = ax_ymin - 2.5*cbar_height
        cbar_left = ax_xmin + 0.5*ax_delta_x - 0.5*cbar_width
        
        cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))        
        cbar = plt.colorbar(im, cax=cax, orientation='horizontal')
        
        # make a "title" (label "m/s" to the right of the colorbar)
        if logscale:
            locator = ticker.LogLocator(subs='all')
            cbar.set_ticks(locator)
            cbar_units = ' ' + texunits[varname]
        elif posdef:
            cbar_units = ' ' + (r'$\times10^{%i}$' %maxabs_exp) + ' ' +\
                texunits[varname]
            cbar.set_ticks([minmax[0], minmax[1]])
            cbar.set_ticklabels(['%1.1f' %minmax[0], '%1.1f' %minmax[1]])
        elif symlog:
            cbar_units = ' ' + texunits[varname]
            cbar.set_ticks([-minmax[1], -linthresh, 0, linthresh,\
                    minmax[1]])
            cbar.set_ticklabels([sci_format(-minmax[1]),\
                    sci_format(-linthresh), '0', sci_format(linthresh),\
                    sci_format(minmax[1])])
#            cax.minorticks_on()
        else:
            cbar_units = ' ' + (r'$\times10^{%i}$' %maxabs_exp) +\
                    ' ' + texunits[varname]
            cbar.set_ticks([minmax[0], 0, minmax[1]])
            cbar.set_ticklabels(['%1.1f' %minmax[0], '0', '%1.1f'\
                    %minmax[1]])
        # Title the colorbar based on the field's units
        fig.text(cbar_left + cbar_width, cbar_bottom + 0.5*cbar_height,\
                 cbar_units, verticalalignment='center', **csfont,\
                 fontsize=cbar_fs) 
    
    # Plot outer boundary
    psivals = np.linspace(0, 2*np.pi, 100)
    xvals, yvals = 2.*np.cos(psivals), np.sin(psivals)
    ax.plot(xvals, yvals, 'k', linewidth=1.3*lw_scaling)
    if figwasNone: # user probably called plot_moll from the python 
        # command line, wanting to view the projection immediately
        plt.show()