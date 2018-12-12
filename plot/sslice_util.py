import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('/altair/loma3853/rayleigh/plot')
from binormalized_cbar import MidpointNormalize
from mpl_toolkits.basemap import Basemap, addcyclic
from matplotlib import colors

def rms(array):
    if np.size(array) == 0:
        return 0
    else:
        return np.sqrt(np.mean(array**2))

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

var_indices = {\
    'vr'    :       1, 
    'vt'    :       2,
    'vp'    :       3,
    'omr'   :       103,
    'omt'   :       104,
    'omp'   :       105 }

labels = {\
    'vr'        :       r'$v_r$',
    'vt'        :       r'$v_\theta$',
    'vp'        :       r'$v_\phi$',
    'vp_prime'  :       r'$v_\phi - \langle v_\phi\rangle$',
    'vl'        :       r'$v_\lambda$',
    'vz'        :       r'$v_z$',
    'omr'       :       r'$\omega_r$',
    'omt'       :       r'$\omega_\theta$',
    'omp'       :       r'$\omega_\phi$',
    'oml'       :       r'$\omega_\lambda$',
    'omz'       :       r'$\omega_z$',
    'vsq'       :       r'$v^2$',
    'vhsq'      :       r'$v_h^2$',
    'vrsq'      :       r'$v_r^2$',
    'omsq'      :       r'$\omega^2$',
    's'         :       r'$S - \langle S \rangle_s$',
    
    'vrs'       :       r'$v_r(S - \langle S\rangle_s)$',
    'vrvp'      :       r'$v_r(v_\phi - \langle v_\phi\rangle)$',
    'vrvt'      :       r'$v_rv_\theta$',
    'vtvp'      :       r'$v_\theta (v_\phi - \langle v_\phi\rangle)$',
    
    'vlvp'      :       r'$v_\lambda (v_\phi - \langle v_\phi\rangle)$',
    'vzvp'      :       r'$v_z (v_\phi - \langle v_\phi\rangle)$',
    'vlvz'      :       r'$v_\lambda v_z$' }

units = {\
    'vr'        :       r'$\frac{\rm{m}}{\rm{s}}$',
    'vt'        :       r'$\frac{\rm{m}}{\rm{s}}$',
    'vp'        :       r'$\frac{\rm{m}}{\rm{s}}$',
    'vp_prime'  :       r'$\frac{\rm{m}}{\rm{s}}$',
    'vl'        :       r'$\frac{\rm{m}}{\rm{s}}$',
    'vz'        :       r'$\frac{\rm{m}}{\rm{s}}$',

    'omr'       :       r'$\frac{\rm{rad}}{\rm{s}}$',
    'omt'       :       r'$\frac{\rm{rad}}{\rm{s}}$',
    'omp'       :       r'$\frac{\rm{rad}}{\rm{s}}$',
    'oml'       :       r'$\frac{\rm{rad}}{\rm{s}}$',
    'omz'       :       r'$\frac{\rm{rad}}{\rm{s}}$',

    'vsq'       :       r'$\frac{\rm{m}^2}{\rm{s}^2}$',
    'vhsq'      :       r'$\frac{\rm{m}^2}{\rm{s}^2}$',
    'vrsq'      :       r'$\frac{\rm{m}^2}{\rm{s}^2}$',
    'omsq'      :       r'$\frac{\rm{rad}^2}{\rm{s}^2}$',
    's'         :       r'$\frac{\rm{erg}}{\rm{K}\ \rm{g}}$',
    
    'vrs'       :       r'$\frac{\rm{erg}}{\rm{K}\ \rm{g}}\frac{\rm{m}}{\rm{s}}$',
    'vrvp'      :       r'$\frac{\rm{m}^2}{\rm{s}^2}$',
    'vrvt'      :       r'$\frac{\rm{m}^2}{\rm{s}^2}$',
    'vtvp'      :       r'$\frac{\rm{m}^2}{\rm{s}^2}$',
    
    'vlvp'      :       r'$\frac{\rm{m}^2}{\rm{s}^2}$',
    'vzvp'      :       r'$\frac{\rm{m}^2}{\rm{s}^2}$',
    'vlvz'      :       r'$\frac{\rm{m}^2}{\rm{s}^2}$' }