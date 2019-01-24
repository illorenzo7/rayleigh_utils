import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['co'])
sys.path.append(os.environ['rapp'])
from binormalized_cbar import MidpointNormalize
from mpl_toolkits.basemap import Basemap, addcyclic
from matplotlib import colors
from varprops import texlabels, texunits, var_indices, var_indices_old
from common import get_widest_range_file
from get_parameter import get_parameter
from rayleigh_diagnostics import ReferenceState

def rms(array):
    if np.size(array) == 0:
        return 0
    else:
        return np.sqrt(np.mean(array**2))

def show_ortho(field, varname, minmax=None, ri=5, ro=6.586, rloc = 6.586, iterr=0):
    # Get latitude and longitude grid from the shape of field
    nphi, nt = np.shape(field)
    dlon = 360.0/nphi
    dlat = 180.0/(nt + 1)
    lons = np.zeros(nphi)
    for i in range(nphi):
        lons[i] = dlon*i-180.0

    lats = np.zeros(nt)
    for i in range(nt):
        lats[i] = (i+1)*dlat - 90    
    
    posdef = False
    if 'sq' in varname:
        posdef = True
        
    if not posdef:
        rms_plus = rms(field[np.where(field > 0)])
        rms_minus = rms(field[np.where(field < 0)])
        if minmax is None:
            my_min, my_max = -3*rms_minus, 3*rms_plus
        else:
            my_min, my_max = minmax
        minexp = int(np.floor(np.log10(np.abs(my_min))))
        maxexp = int(np.floor(np.log10(np.abs(my_max))))
        maxabs_exp = max((minexp, maxexp))

        field /= 10**maxabs_exp
        my_min /= 10**maxabs_exp
        my_max /= 10**maxabs_exp
        
    # Set up the lat/lon grid for the Basemap projections
    xpixels = 512 
    ypixels = 512

    # Set up the cyclic grid to avoid "cuts" of whitespace between 
    # 0 and 360 degrees
    sslice_cycl, lons_cycl = addcyclic(field, lons)

    # Convert to 2-D grids
    llons_cycl, llats_cycl = np.meshgrid(lons_cycl, lats)
    
    # Create the plot using subplot axes
    # Offset axes slightly (at the end) to deal with annoying white space cutoff
    offset_fraction = 0.01
    fig_width_inches = 8.
    
    # General parameters for main axis/color bar
    margin_bottom_inches = 1.
    margin_top_inches = 1/2
    margin_inches = 1/8

    subplot_width_inches = fig_width_inches - 2*margin_inches
    subplot_height_inches = subplot_width_inches
    fig_height_inches = margin_bottom_inches + subplot_height_inches + margin_top_inches
    fig_aspect = fig_height_inches/fig_width_inches
    
    # "Non-dimensional" figure parameters
    margin_bottom = margin_bottom_inches/fig_height_inches
    margin_top = margin_top_inches/fig_height_inches
    subplot_width = subplot_width_inches/fig_width_inches
    subplot_height = subplot_height_inches/fig_height_inches
    margin_x = margin_inches/fig_width_inches
    margin_y = margin_inches/fig_height_inches
    
    # Main axis plot: orthographic projection
    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax = fig.add_axes([margin_x, margin_bottom, subplot_width, subplot_height])
#    plt.sca(ax)
    m = Basemap(projection='ortho', lon_0=0, lat_0=20, resolution=None, rsphere=1.)
    field = np.transpose(field)
    topodat, x, y = m.transform_scalar(field, lons, lats, xpixels, ypixels, returnxy=True, masked=True, order=1)
    
    shrink_distance = 1. - rloc/ro
    shrink_factor = rloc/ro
    d = 1. - ri/ro
    depth = (1. - rloc/ro)/d    
    x = x*shrink_factor + shrink_distance
    y = y*shrink_factor + shrink_distance
    
    # Saturation levels for a positive-definite quantity (like vsq)
    if posdef:
        logfield= np.log(field)
        medlog = np.median(logfield)
        shiftlog = logfield - medlog
        minexp = medlog - 7*np.std(shiftlog[np.where(shiftlog < 0)].flatten())
        maxexp = medlog + 7*np.std(shiftlog[np.where(shiftlog > 0)].flatten())
        my_min, my_max = np.exp(minexp), np.exp(maxexp)
        
    #View the data
    if not posdef:
        m.pcolormesh(x, y, topodat, cmap=plt.cm.RdYlBu_r,\
                     norm=MidpointNormalize(0), vmin=my_min, vmax=my_max)
    else: 
        m.pcolormesh(x,y, topodat, cmap='Greys',\
            norm=colors.LogNorm(vmin=my_min, vmax=my_max))
        
    # draw parallels and meridians, evenly spaced by 30 degrees
    parallels = np.arange(-60, 90, 30.)
    meridians = np.arange(0., 360., 30.)
    
    for meridian in meridians:
        linex, liney = m(np.ones_like(lats)*meridian, lats)

        # must cut off all points on the far side of the sphere!
        linex = linex[np.where(np.abs(linex) < 1.e20)]
        liney = liney[np.where(np.abs(liney) < 1.e20)]

        linex = linex*shrink_factor + shrink_distance
        liney = liney*shrink_factor + shrink_distance

        m.plot(linex, liney, color='k', linewidth=.5)
        
    for parallel in parallels:
        linex, liney = m(lons, np.ones_like(lons)*parallel)

        # must cut off all points on the far side of the sphere!
        linex = linex[np.where(np.abs(linex) < 1.e20)]
        liney = liney[np.where(np.abs(liney) < 1.e20)]

        linex = linex*shrink_factor + shrink_distance
        liney = liney*shrink_factor + shrink_distance

        m.plot(linex, liney, color='k', linewidth=.5)        
    
    # Now adjust axes slightly to get rid of whitespace
    xmin, xmax = ax.get_xlim(); delta_x = xmax - xmin
    ax.set_xlim((xmin - offset_fraction*delta_x, xmax + offset_fraction*delta_x))
    ymin, ymax = ax.get_ylim(); delta_y = ymax - ymin
    ax.set_ylim((ymin - offset_fraction*delta_y, ymax + offset_fraction*delta_y))
    
    # Set up color bar
    cbar_aspect = 1./20.
    cbarregion_center_x = margin_x + 0.5*subplot_width
    cbarregion_center_y = 0.5*margin_bottom
    cbar_width = 0.5*subplot_width
    cbar_height = cbar_aspect*cbar_width/fig_aspect
    cbar_left = cbarregion_center_x - 0.5*cbar_width
    cbar_bottom = margin_bottom - 1.2*cbar_height
    cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))
    
    cbar = plt.colorbar(cax=cax, orientation='horizontal')
    # make a "title" (label "m/s" to the right of the colorbar)
    if not posdef:
        cbar_units = ' ' + (r'$\times10^{%i}$' %maxabs_exp) + ' ' + texunits[varname]
        cbar.set_ticks([my_min, 0, my_max])
        cbar.set_ticklabels(['%1.2f' %my_min, '0', '%1.2f' %my_max])
    else:
        cbar_units = ' ' + texunits[varname]
    varlabel = texlabels[varname]

    title = varlabel + '\t\t' + ('depth = %1.2f' %depth) + '\t\t' + ('iter = %08i' %iterr)
    fig.text(cbar_left + cbar_width, cbar_bottom + 0.5*cbar_height, cbar_units,\
                 verticalalignment='center')
    fig.text(margin_x + 0.5*subplot_width, 1. - 0.6*margin_top, title,\
             verticalalignment='center', horizontalalignment='center', fontsize=14)    

    plt.show()   
    
def plot_ortho(fig, ax, field, lats, lons, r_loc, ri, ro,\
               cax=None, posdef=False, cbar_units = 'm/s', minmax='default',
               return_minmax = False):
    plt.sca(ax)
    
    use_sci = False
    if not posdef:
        rms_plus = rms(field[np.where(field > 0)])
        rms_minus = rms(field[np.where(field < 0)])
        my_min, my_max = -3*rms_minus, 3*rms_plus
        minexp = int(np.floor(np.log10(np.abs(my_min))))
        maxexp = int(np.floor(np.log10(np.abs(my_min))))
        maxabs_exp = max((minexp, maxexp))

        if (maxabs_exp > 2):
            use_sci = True
            field /= 10**maxabs_exp
            my_min /= 10**maxabs_exp
            my_max /= 10**maxabs_exp
        
    # Set up the lat/lon grid for the Basemap projections
    xpixels = 1024  
    ypixels = 1024
    
    # Compute the location of the tangent cylinder at this depth
    tangent_theta = np.arcsin(ri/r_loc)
    tangent_theta_deg = 180./np.pi*tangent_theta
    tangent_lat = 90. - tangent_theta_deg

    # Set up the cyclic grid for Mollweide projection
    sslice_cycl, lons_cycl = addcyclic(field, lons)

    # Convert to 2-D grids
    llons_cycl, llats_cycl = np.meshgrid(lons_cycl, lats)
    
    # Make the Mollweide plot
    #fig = plt.figure(figsize=(6,3), dpi=300)
    m = Basemap(projection='ortho', lon_0=0, lat_0=20, resolution=None, rsphere=ro)

    topodat,x,y =\
    m.transform_scalar(field, lons, lats, xpixels, ypixels, \
                       returnxy=True, masked=True, order=1)

    shrink_distance = ro - r_loc
    shrink_factor = ro/r_loc
    x = x/shrink_factor + shrink_distance
    y = y/shrink_factor + shrink_distance 

    # Saturation levels for a positive-definite quantity (like vsq)
    if posdef:
        logfield= np.log(field)
        medlog = np.median(logfield)
        shiftlog = logfield - medlog
        minexp = medlog - 7*np.std(shiftlog[np.where(shiftlog < 0)].flatten())
        maxexp = medlog + 7*np.std(shiftlog[np.where(shiftlog > 0)].flatten())
        my_min, my_max = np.exp(minexp), np.exp(maxexp)
        
    if not minmax == 'default':
        my_min, my_max = minmax
        
    #View the data
    if not posdef:
        m.pcolormesh(x, y, topodat, cmap=plt.cm.RdYlBu_r,\
                     norm=MidpointNormalize(0), vmin=my_min, vmax=my_max)
    else: 
        m.pcolormesh(x,y, topodat, cmap='Greens',\
            norm=colors.LogNorm(vmin=my_min, vmax=my_max))
    # draw parallels and meridians. Draw two parallels at the tangent cylinder
#    m.drawparallels((-tangent_lat, tangent_lat))
#    m.drawmeridians(np.arange(0.,420.,60.))
    
    # convert desires lons/lats to draw into map projection coordinates
    meridians = np.arange(0., 360., 60.)
    
    for meridian in meridians:
        linex, liney = m(np.ones_like(lats)*meridian, lats)
    
        # must cut off all points on the far side of the sphere!
        linex = linex[np.where(np.abs(linex) < 1.e20)]
        liney = liney[np.where(np.abs(liney) < 1.e20)]
    
        linex = linex/shrink_factor + shrink_distance
        liney = liney/shrink_factor + shrink_distance
    
        m.plot(linex, liney,color='k',linewidth=.75)
    
    # now plot where the tangent cylinder intersects the shell slice:
    for tlat in [-tangent_lat, tangent_lat]:
        linex, liney = m(lons, np.ones_like(lons)*tlat)
        # must cut off all points on the far side of the sphere!
        linex = linex[np.where(np.abs(linex) < 1.e20)]
        liney = liney[np.where(np.abs(liney) < 1.e20)]
    
        linex = linex/shrink_factor + shrink_distance
        liney = liney/shrink_factor + shrink_distance
    
        m.plot(linex, liney,'k',linewidth=.75)
        
    if not cax is None:
        cbar = plt.colorbar(cax=cax, orientation='horizontal')
        # make a "title" (label "m/s" to the right of the colorbar)
        xmin, xmax = cbar.ax.get_xlim(); delta_x = xmax - xmin
        ymin, ymax = cbar.ax.get_ylim(); delta_y = ymax - ymin
        cbar_label = cbar_units
        
        if use_sci:
            cbar_label = (r'$\times10^{%i}\ $' %maxabs_exp) + cbar_units
        
        cbar.ax.text(xmax + 0.15*delta_x, ymin + delta_y/2, cbar_label,\
                     verticalalignment='center', fontsize=12)
        if not posdef:
            cbar.set_ticks([my_min, 0, my_max])
        if use_sci:
            cbar.set_ticklabels(['%1.1f' %my_min, '0', '%1.1f' %my_max])
        else:
            if not posdef:
                cbar.set_ticklabels(['%2.0f' %my_min, '0', '%2.0f' %my_max])
    if return_minmax:
        return my_min, my_max

    
def cbar_axes_from_area(cbar_area_left, cbar_area_bottom, \
                        cbar_area_width, cbar_area_height, fig_aspect, aspect=1/20):
    cbar_area_centerx = cbar_area_left + cbar_area_width/2
    cbar_area_centery = cbar_area_bottom + cbar_area_height/2
    cbar_height = cbar_area_height/8
    cbar_width = cbar_height/aspect*fig_aspect
    cbar_left = cbar_area_centerx - cbar_width/2
    cbar_bottom = cbar_area_centery + cbar_area_height/4 - cbar_height/2
    return (cbar_left, cbar_bottom, cbar_width, cbar_height)