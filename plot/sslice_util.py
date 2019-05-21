import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib import ticker
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['co'])
sys.path.append(os.environ['rapp'])
from binormalized_cbar import MidpointNormalize
from cartopy import crs
from matplotlib import colors
from varprops import texlabels, texunits, var_indices, var_indices_old
from common import get_file_lists, get_satvals,\
        saturate_array, rsun
from get_sslice import get_sslice
from rayleigh_diagnostics import Shell_Slices
from get_parameter import get_parameter

def deal_with_nans(x, y):
    """handles NaN in x and y by making quadrilateral corners on the "other
    side" of the projection overlap with adjacent quadrilateral corners--so
    that pcolormesh plots nothing there
    """
    x_dealt, y_dealt = np.copy(x), np.copy(y)

    mask = np.isfinite(x) & np.isfinite(y)
    
    found_min_iy = False
    found_max_iy = False
    
    nx, ny = np.shape(mask)
    
    for iy in range(ny):
        good = mask[:, iy].nonzero()[0]
        # Break "good" into blocks that are individually contiguous       
        if good.size == 0:
            continue
        elif not found_min_iy:
            min_iy = iy
            found_min_iy = True
        else:
            max_iy = iy
            found_max_iy = True
        dgood = np.diff(good)
        if True in (dgood > 1): # The viewable longitudes cross over 180
            # so the bad values occur as a chunk in the middle of the array
            len_notgood = dgood[np.where(dgood > 1)][0] - 1
            start_notgood = np.where(dgood > 1)[0][0] + 1
            x_dealt[start_notgood:start_notgood + len_notgood, iy] =\
                x_dealt[start_notgood - 1, iy]
            y_dealt[start_notgood:start_notgood + len_notgood, iy] =\
                y_dealt[start_notgood - 1, iy]                
        else:
            x_dealt[good[-1]:, iy] = x_dealt[good[-1], iy]
            y_dealt[good[-1]:, iy] = y_dealt[good[-1], iy]

            x_dealt[:good[0], iy] = x_dealt[good[0], iy]
            y_dealt[:good[0], iy] = y_dealt[good[0], iy]
    
    if found_min_iy:
        for ix in range(nx):
            x_dealt[ix, :min_iy] = x_dealt[ix, min_iy]
            y_dealt[ix, :min_iy] = y_dealt[ix, min_iy]

    if found_max_iy:
        for ix in range(nx):
            x_dealt[ix, max_iy:] = x_dealt[ix, max_iy]
            y_dealt[ix, max_iy:] = y_dealt[ix, max_iy]
    
    radius = np.max(np.sqrt(x_dealt**2 + y_dealt**2))
    return x_dealt/radius, y_dealt/radius

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
    
#def get_ortho_parallels(lat, pc, ortho):
    
def plot_ortho(fig, ax, a, dirname, varname, ir=0, minmax=None,\
               clon=0, clat=20, plot_title=True):
    
    fname = str(a.iters[0]).zfill(8)
    vals = get_sslice(a, varname, dirname=dirname)
    field = vals[:, :, ir]
    
    # Get geometric parameters
    ri, ro = rbounds(dirname)
    rloc = a.radius[ir]
    
    # Get latitude and longitude grid from the shape of field
    nphi, ntheta = np.shape(field)
    dlon = 360.0/nphi
    dlat = 180.0/ntheta
    lons = dlon*np.arange(nphi) - 180.0
    lats = dlat*np.arange(ntheta) - 90.0
    # Make 2d (meshed) grid of longitude/latitude
    llon, llat = np.meshgrid(lons, lats, indexing='ij')
    
    # Essence of cartopy is transformations:
    # This one will be from PlateCarree (simple longitude/latitude equally 
    # spaced) to orthographic
    pc = crs.PlateCarree()
    ortho = crs.Orthographic(central_latitude=clat, central_longitude=clon)
    
    # Determinate if "field" is positive-definite (inherently)
    posdef = False
    if 'sq' in varname:
        posdef = True
    
    # Get saturation values to be used if minmax is not specified
    # Divide out the exponent to use scientific notation (but keep 
    # track of it!)
    if minmax is None:
        my_min, my_max = get_satvals(field, posdef=posdef)
        if not posdef: # normalize values to plot in scientific notation
            minexp = int(np.floor(np.log10(np.abs(my_min))))
            maxexp = int(np.floor(np.log10(np.abs(my_max))))
            maxabs_exp = max((minexp, maxexp))
    
            field /= 10**maxabs_exp
            my_min /= 10**maxabs_exp
            my_max /= 10**maxabs_exp    

    else: # minmax IS none
            my_min, my_max = minmax
            minexp = int(np.floor(np.log10(np.abs(my_min))))
            maxexp = int(np.floor(np.log10(np.abs(my_max))))
            maxabs_exp = max((minexp, maxexp))
    
            field /= 10**maxabs_exp
            my_min /= 10**maxabs_exp
            my_max /= 10**maxabs_exp              
            
    # Get the orthographic projection coordinates of llon/llat by converting
    # between PlateCarree (lat/lon) and orthographic
    x_withnan = np.zeros((nphi, ntheta)) 
    y_withnan = np.zeros((nphi, ntheta)) 
    for iphi in range(nphi):
        for itheta in range(ntheta):
            lon_loc, lat_loc = llon[iphi, itheta], llat[iphi, itheta]
            x_withnan[iphi, itheta], y_withnan[iphi, itheta] =\
                ortho.transform_point(lon_loc, lat_loc, pc)
    # Must deal with "nan" coordinates since part of sphere is invisible
    x_unshrunk, y_unshrunk = deal_with_nans(x_withnan, y_withnan)
    
    # May possibly need to shrink coordinates to give depth perception
    shrink_factor = rloc/ro
    d = 1. - ri/ro
    depth = (1. - rloc/ro)/d    
    x = x_unshrunk*shrink_factor #+ shrink_distance
    y = y_unshrunk*shrink_factor #+ shrink_distance
    

    ax.set_xlim((-1.01, 1.01)) # deal with annoying whitespace cutoff issue
    ax.set_ylim((-1.01, 1.01))
    ax.axis('off') # get rid of x/y axis coordinates
            
    # Make the orthographic projection
    if not posdef:
        saturate_array(field, my_min, my_max)
#        im = ax.pcolormesh(x, y, field, cmap=plt.cm.RdYlBu_r,\
#                     norm=MidpointNormalize(0), vmin=my_min, vmax=my_max)
        im = ax.contourf(x, y, field, cmap=plt.cm.RdYlBu_r,\
                levels=np.linspace(my_min, my_max, 50),\
                norm=MidpointNormalize(0))
    else: 
#        im = ax.pcolormesh(x, y, field, cmap='Greys',\
#            norm=colors.LogNorm(vmin=my_min, vmax=my_max))
         im = ax.contourf(x, y, field, cmap='Greys',\
            norm=colors.LogNorm(vmin=my_min, vmax=my_max),\
            levels=np.logspace(minexp, maxexp, 50, base=np.exp(1.)))
       
    # Draw parallels and meridians, evenly spaced by 30 degrees
    default_lw = 0.5 # default linewidth bit thinner
    parallels = np.arange(-60, 90, 30.)
    meridians = np.arange(0., 360., 30.)
    
    npoints = 100
    for meridian in meridians:
        if meridian == 0.: 
            lw = 1.3 # make central longitude thicker
        else:
            lw = default_lw
            
        lonvals = meridian + np.zeros(npoints)
        latvals = np.linspace(-90, 90, npoints)
        linex_withnans, liney_withnans = np.zeros(npoints), np.zeros(npoints)
        for i in range(npoints):
            linex_withnans[i], liney_withnans[i] =\
                ortho.transform_point(lonvals[i], latvals[i], pc)
        linex, liney = deal_with_nans(linex_withnans.reshape((npoints, 1)),\
                                      liney_withnans.reshape((npoints, 1)))
        linex = linex[:, 0]*shrink_factor
        liney = liney[:, 0]*shrink_factor
        ax.plot(linex, liney, 'k', linewidth=lw)
    
    for parallel in parallels:
        if parallel == 0.: 
            lw = 1.3 # make equator thicker
        else:
            lw = default_lw        
        lonvals = np.linspace(-180, 180, npoints)
        latvals = np.zeros(npoints) + parallel
        linex_withnans, liney_withnans = np.zeros(npoints), np.zeros(npoints)
        for i in range(npoints):
            linex_withnans[i], liney_withnans[i] =\
                ortho.transform_point(lonvals[i], latvals[i], pc)
        linex, liney = deal_with_nans(linex_withnans.reshape((npoints, 1)),\
                                      liney_withnans.reshape((npoints, 1)))
        linex = linex[:, 0]*shrink_factor
        liney = liney[:, 0]*shrink_factor
        ax.plot(linex, liney, 'k', linewidth=lw)   

    # Set up color bar
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x
#    ax_center_y = ax_ymin + 0.5*ax_delta_y   
    
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
    # make a "title" (label "m/s" to the right of the colorbar)
    if not posdef:
        cbar_units = ' ' + (r'$\times10^{%i}$' %maxabs_exp) + ' ' + texunits[varname]
        cbar.set_ticks([my_min, 0, my_max])
        cbar.set_ticklabels(['%1.1f' %my_min, '0', '%1.1f' %my_max])
    else:
        locator = ticker.LogLocator(base=10)
        cbar.set_ticks(locator)
        cbar_units = ' ' + texunits[varname]
        
    varlabel = texlabels[varname]

    title = varlabel + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %(rloc/rsun)) +\
            '     ' + ('iter = ' + fname)
    fig.text(cbar_left + cbar_width, cbar_bottom + 0.5*cbar_height,\
             cbar_units, verticalalignment='center', **csfont)
    if plot_title:
        fig.text(ax_center_x, ax_ymax + 0.02*ax_delta_y, title,\
             verticalalignment='bottom', horizontalalignment='center',\
             fontsize=14, **csfont)   
    
    # Plot outer boundary
    psivals = np.linspace(0, 2*np.pi, 100)
    xvals, yvals = np.cos(psivals), np.sin(psivals)
    ax.plot(xvals, yvals, 'k')

def plot_moll(fig, ax, a, dirname, varname, ir=0, minmax=None,\
               clon=0, plot_title=True):
    
    fname = str(a.iters[0]).zfill(8)
    vals = get_sslice(a, varname, dirname=dirname)
    field = vals[:, :, ir]
    
    # Get geometric parameters
    ri, ro = rbounds(dirname)
    rloc = a.radius[ir]
    d = 1. - ri/ro
    depth = (1. - rloc/ro)/d    
    
    # Get latitude and longitude grid from the shape of field
    nphi, ntheta = np.shape(field)
    dlon = 360.0/nphi
    dlat = 180.0/ntheta
    lons = dlon*np.arange(nphi) - 180.0
    lats = dlat*np.arange(ntheta) - 90.0
    # Make 2d (meshed) grid of longitude/latitude
    llon, llat = np.meshgrid(lons, lats, indexing='ij')
    
    # Essence of cartopy is transformations:
    # This one will be from PlateCarree (simple longitude/latitude equally 
    # spaced) to Mollweide
    pc = crs.PlateCarree()
    moll = crs.Mollweide()
    
    # Determinate if "field" is positive-definite (inherently)
    posdef = False
    if 'sq' in varname:
        posdef = True
    
    # Get saturation values to be used if minmax is not specified
    # Divide out the exponent to use scientific notation (but keep 
    # track of it!)
    if minmax is None:
        my_min, my_max = get_satvals(field, posdef=posdef)
    else: # minmax WAS specified by user
        my_min, my_max = minmax

    if not posdef: # normalize values to plot in scientific notation
        minexp = int(np.floor(np.log10(np.abs(my_min))))
        maxexp = int(np.floor(np.log10(np.abs(my_max))))
        maxabs_exp = max((minexp, maxexp))

        field /= 10**maxabs_exp
        my_min /= 10**maxabs_exp
        my_max /= 10**maxabs_exp    
            
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

    ax.set_xlim((-2.02, 2.02)) # deal with annoying whitespace cutoff issue
    ax.set_ylim((-1.01, 1.01))
    ax.axis('off') # get rid of x/y axis coordinates
            
    # Make the Mollweide projection
    if not posdef:
        saturate_array(field, my_min, my_max)
        im = ax.contourf(x, y, field, cmap=plt.cm.RdYlBu_r,\
                levels=np.linspace(my_min, my_max, 50),\
                norm=MidpointNormalize(0))
    else: 
         im = ax.contourf(x, y, field, cmap='Greys',\
            norm=colors.LogNorm(vmin=my_min, vmax=my_max),\
            levels=np.logspace(minexp, maxexp, 50, base=np.exp(1.)))
       
    # Draw parallels and meridians, evenly spaced by 30 degrees
    default_lw = 0.5 # default linewidth bit thinner
    parallels = np.arange(-60, 90, 30.)
    meridians = np.arange(0., 360., 30.)
    
    npoints = 100
    for meridian in meridians:
        if meridian == 0.: 
            lw = 1.3 # make central longitude thicker
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
            lw = 1.3 # make equator thicker
        else:
            lw = default_lw        
        lonvals = np.linspace(-180, 180, npoints)
        latvals = np.zeros(npoints) + parallel
        linepoints = moll.transform_points(pc, lonvals, latvals)
        linex, liney = linepoints[:, 0], linepoints[:, 1]
        linex /= radius
        liney /= radius
        
        ax.plot(linex, liney, 'k', linewidth=lw)

    # Set up color bar
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x
#    ax_center_y = ax_ymin + 0.5*ax_delta_y   
    
    cbar_aspect = 1./20.
    fig_aspect = ax_delta_x/ax_delta_y*0.5 
    # assumes subplot aspect ratio is 0.5
    cbar_width = 0.25*ax_delta_x # make cbar one-quarter
                        # as long as plot is wide
    cbar_height = cbar_width*cbar_aspect/fig_aspect
    cbar_bottom = ax_ymin - 2.5*cbar_height
    cbar_left = ax_xmin + 0.5*ax_delta_x - 0.5*cbar_width
    
    cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))
    
    cbar = plt.colorbar(im, cax=cax, orientation='horizontal')
    # make a "title" (label "m/s" to the right of the colorbar)
    if not posdef:
        cbar_units = ' ' + (r'$\times10^{%i}$' %maxabs_exp) + ' ' + texunits[varname]
        cbar.set_ticks([my_min, 0, my_max])
        cbar.set_ticklabels(['%1.1f' %my_min, '0', '%1.1f' %my_max])
    else:
        locator = ticker.LogLocator(base=10)
        cbar.set_ticks(locator)
        cbar_units = ' ' + texunits[varname]
        
    varlabel = texlabels[varname]

    title = varlabel + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %(rloc/rsun)) +\
            '     ' + ('iter = ' + fname)
    fig.text(cbar_left + cbar_width, cbar_bottom + 0.5*cbar_height,\
             cbar_units, verticalalignment='center', **csfont)
    if plot_title:
        fig.text(ax_center_x, ax_ymax + 0.02*ax_delta_y, title,\
             verticalalignment='bottom', horizontalalignment='center',\
             fontsize=14, **csfont)   
    
    # Plot outer boundary
    psivals = np.linspace(0, 2*np.pi, 100)
    xvals, yvals = 2.*np.cos(psivals), np.sin(psivals)
    ax.plot(xvals, yvals, 'k')
