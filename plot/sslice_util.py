import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['co'])
sys.path.append(os.environ['rapp'])
from binormalized_cbar import MidpointNormalize
from cartopy import crs
from matplotlib import colors
from varprops import texlabels, texunits, var_indices, var_indices_old
from common import get_file_lists
from get_sslice import get_sslice
from rayleigh_diagnostics import Shell_Slices
from get_parameter import get_parameter

def rms(array):
    if np.size(array) == 0:
        return 0
    else:
        return np.sqrt(np.mean(array**2))

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
        if good.size == 0:
            continue
        elif not found_min_iy:
            min_iy = iy
            found_min_iy = True
        else:
            max_iy = iy
            found_max_iy = True

        x_dealt[good[-1]:, iy] = x_dealt[good[-1], iy]
        y_dealt[good[-1]:, iy] = y_dealt[good[-1], iy]

        x_dealt[:good[0], iy] = x_dealt[good[0], iy]
        y_dealt[:good[0], iy] = y_dealt[good[0], iy]
    
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
    rmin = get_parameter(dirname, 'rmin')
    rmax = get_parameter(dirname, 'rmax')
    if rmin == 100 or rmax == 100: # get_parameter must have failed
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
    
def plot_ortho(fig, ax, dirname, varname, idepth=0, minmax=None, iiter=-1):
    
    # Get shell slice desired by caller
    radatadir = dirname + '/Shell_Slices/'
    file_list, int_file_list, nfiles = get_file_lists(radatadir)
    fname = file_list[iiter]
    a = Shell_Slices(radatadir + fname, '')
    vals = get_sslice(a, varname, dirname=dirname)
    field = vals[:, :, idepth]
    
    # Get geometric parameters
    ri, ro = rbounds(dirname)
    rloc = a.radius[idepth]
    
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
    ortho = crs.Orthographic()
    
    # Determinate if "field" is positive-definite (inherently)
    posdef = False
    if 'sq' in varname:
        posdef = True
    
    # Get saturation values to be used if minmax is not specified
    # Divide out the exponent to use scientific notation (but keep 
    # track of it!)
    if minmax is None:
        if not posdef:
            rms_plus = rms(field[np.where(field > 0)])
            rms_minus = rms(field[np.where(field < 0)])
            my_min, my_max = -3*rms_minus, 3*rms_plus
            
            minexp = int(np.floor(np.log10(np.abs(my_min))))
            maxexp = int(np.floor(np.log10(np.abs(my_max))))
            maxabs_exp = max((minexp, maxexp))
    
            field /= 10**maxabs_exp
            my_min /= 10**maxabs_exp
            my_max /= 10**maxabs_exp    
            
        # Saturation levels for a positive-definite quantity (like vsq)
        else:
            logfield= np.log(field)
            medlog = np.median(logfield)
            shiftlog = logfield - medlog
            minexp = medlog - 7*np.std(shiftlog[np.where(shiftlog < 0)].flatten())
            maxexp = medlog + 7*np.std(shiftlog[np.where(shiftlog > 0)].flatten())
            my_min, my_max = np.exp(minexp), np.exp(maxexp)
    else: # minmax IS none
            my_min, my_max = minmax
            
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
    shrink_distance = 1. - rloc/ro
    shrink_factor = rloc/ro
    d = 1. - ri/ro
    depth = (1. - rloc/ro)/d    
    x = x_unshrunk*shrink_factor #+ shrink_distance
    y = y_unshrunk*shrink_factor #+ shrink_distance
    

    ax.set_xlim((-1.01, 1.01)) # deal with annoying whitespace cutoff issue
    ax.set_ylim((-1.01, 1.01))
    plt.axis('off') # get rid of x/y axis coordinates
            
    # Make the orthographic projection
    if not posdef:
        im = ax.pcolormesh(x, y, field, cmap=plt.cm.RdYlBu_r,\
                     norm=MidpointNormalize(0), vmin=my_min, vmax=my_max)
    else: 
        im = ax.pcolormesh(x,y, field, cmap='Greys',\
            norm=colors.LogNorm(vmin=my_min, vmax=my_max))
        
    # draw parallels and meridians, evenly spaced by 30 degrees
    parallels = np.arange(-60, 90, 30.)
    meridians = np.arange(0., 360., 30.)
    
#    for meridian in meridians:
#        linex, liney = m(np.ones_like(lats)*meridian, lats)
#
#        # must cut off all points on the far side of the sphere!
#        linex = linex[np.where(np.abs(linex) < 1.e20)]
#        liney = liney[np.where(np.abs(liney) < 1.e20)]
#
#        linex = linex*shrink_factor + shrink_distance
#        liney = liney*shrink_factor + shrink_distance
#
#        m.plot(linex, liney, color='k', linewidth=.5)
#        
#    for parallel in parallels:
#        linex, liney = m(lons, np.ones_like(lons)*parallel)
#
#        # must cut off all points on the far side of the sphere!
#        linex = linex[np.where(np.abs(linex) < 1.e20)]
#        liney = liney[np.where(np.abs(liney) < 1.e20)]
#
#        linex = linex*shrink_factor + shrink_distance
#        liney = liney*shrink_factor + shrink_distance
#
#        m.plot(linex, liney, color='k', linewidth=.5)        
    
    # Now adjust axes slightly to get rid of whitespace
    
#    xmin, xmax = ax.get_xlim(); delta_x = xmax - xmin
#    ax.set_xlim((xmin - offset_fraction*delta_x, xmax + offset_fraction*delta_x))
#    ymin, ymax = ax.get_ylim(); delta_y = ymax - ymin
#    ax.set_ylim((ymin - offset_fraction*delta_y, ymax + offset_fraction*delta_y))

   
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
    cbar_bottom = ax_ymin - 1.5*cbar_height
    cbar_left = ax_xmin + 0.5*ax_delta_x - 0.5*cbar_width
    
#    cbarregion_center_x = margin_x + 0.5*subplot_width
#    cbarregion_center_y = 0.5*margin_bottom
#    cbar_width = 0.5*subplot_width
#    cbar_height = cbar_aspect*cbar_width/fig_aspect
#    cbar_left = cbarregion_center_x - 0.5*cbar_width
#    cbar_bottom = margin_bottom - 1.2*cbar_height
    cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))
    
    cbar = plt.colorbar(im, cax=cax, orientation='horizontal')
    # make a "title" (label "m/s" to the right of the colorbar)
    if not posdef:
        cbar_units = ' ' + (r'$\times10^{%i}$' %maxabs_exp) + ' ' + texunits[varname]
        cbar.set_ticks([my_min, 0, my_max])
        cbar.set_ticklabels(['%1.2f' %my_min, '0', '%1.2f' %my_max])
    else:
        cbar_units = ' ' + texunits[varname]
        
    varlabel = texlabels[varname]

    title = varlabel + '\t\t' + ('depth = %1.2f' %depth) + '\t\t' +\
        ('iter = ' + fname)
    fig.text(cbar_left + cbar_width, cbar_bottom + 0.5*cbar_height,\
             cbar_units, verticalalignment='center')
    fig.text(ax_center_x, ax_ymax + 0.02*ax_delta_y, title,\
             verticalalignment='bottom', horizontalalignment='center',\
             fontsize=14)   
    
    # Plot outer boundary
    psivals = np.linspace(0, 2*np.pi, 100)
    xvals, yvals = np.cos(psivals), np.sin(psivals)
    ax.plot(xvals, yvals, 'k')
#    plt.sca(ax)
    plt.show()   

    
def cbar_axes_from_area(cbar_area_left, cbar_area_bottom, \
                        cbar_area_width, cbar_area_height, fig_aspect, aspect=1/20):
    cbar_area_centerx = cbar_area_left + cbar_area_width/2
    cbar_area_centery = cbar_area_bottom + cbar_area_height/2
    cbar_height = cbar_area_height/8
    cbar_width = cbar_height/aspect*fig_aspect
    cbar_left = cbar_area_centerx - cbar_width/2
    cbar_bottom = cbar_area_centery + cbar_area_height/4 - cbar_height/2
    return (cbar_left, cbar_bottom, cbar_width, cbar_height)