# Author: Loren Matilsky
# Created: 03/05/2018
####################################################################################################
#
#  Shell-Slice (Shell_Slices) plotting example
#  - Reads in a single Shell_Slice file.
#  - Plots velocity components in both spherical and cylindrical 
#       coordinates (5 components total), vorticity (5 components, if they 
#       are available), enstrophy, velocity-squared, "horizontal" velocity-
#       squared, "vertical" (vr) velocity-squared, and entropy at all 
#       depths for a given slice. 
#  - In total, script attempts to plot FIFTEEN (15) fluid variables.
# 
#  - Saves high-res .png files (dpi=300) in panel-form (1 panel for each
#       depth), and separately depth-by-depth in the folder 
#       plots/sslices/varname
#  
#  - Saves in both Mollweide and Orthographic projections
#
#  - This example routine makes use of the ShellSlice
#       data structure associated with the Shell_Slices output.
#       Upon initializing a ShellSlice object, the 
#       object will contain the following attributes:
#    ----------------------------------
#    self.niter                 : number of time steps
#    self.nq                    : number of diagnostic quantities output
#    self.nr                    : number of shell slices output
#    self.ntheta                : number of theta points
#    self.nphi                  : number of phi points
#    self.qv[0:nq-1]            : quantity codes for the diagnostics output
#    self.radius[0:nr-1]        : radii of the shell slices output
#    self.inds[0:nr-1]          : radial indices of the shell slices output
#    self.costheta[0:ntheta-1]  : cos(theta grid)
#    self.sintheta[0:ntheta-1]  : sin(theta grid)
#    self.vals[0:nphi-1,0:ntheta-1,0:nr-1,0:nq-1,0:niter-1] 
#                               : The shell slices 
#    self.iters[0:niter-1]      : The time step numbers stored in this 
#                                   output file
#    self.time[0:niter-1]       : The simulation time corresponding to 
#                                   each time step
#    self.version               : The version code for this particular 
#                                   output (internal use)
#    self.lut                   : Lookup table for the different 
#                                   diagnostics output
#    -------------------------------------
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
#from diagnostic_reading import ShellSlice
from rayleigh_diagnostics import Shell_Slices
from matplotlib import ticker
from matplotlib import colors
from mpl_toolkits.basemap import Basemap, addcyclic
from common import get_widest_range_file

# Get the run directory on which to perform the analysis
#dirname = sys.argv[1] # changed for Brad
fname = sys.argv[1]
dirname = 'sim_' + fname
os.makedirs(dirname + '/Shell_Slices')
os.rename(fname, dirname + '/Shell_Slices/' + fname)

# Directory with data
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'

# Read command-line arguments (CLAs)
number_from_end = 1 # By default, plot the last slice
args = sys.argv[2:]
nargs = len(args)
oldra = False
for i in range(nargs):
    arg = args[i]
    if (arg == '-n'):
        number_from_end = int(args[i+1])
    elif (arg == '-plotdir'):
        plotdir = dirname + '/' + args[i+1] + '/'
    elif (arg == '-old'):
        oldra = True

if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Formatting function for colorbar
def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

# Get grid info
#rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
#nt, nr = len(tt), len(rr)
#nphi = 2*nt
#cost_nd = cost.reshape((1, nt))
#sint_nd = sint.reshape((1, nt))

#rsun = 6.96e10
#rr_n = rr/rsun # rr "normalized" by the solar radius

# Get mean entropy
#savg_file = get_widest_range_file(datadir, 's_spherical_mean')
#s_spherical_mean = np.load(datadir + savg_file)

# Get mean velocities
#vavg_file = get_widest_range_file(datadir, 'vavg')
#vavg = np.load(datadir + vavg_file)
#vr_av, vt_av, vp_av = vavg


# Get a list of all the shellslices
radatadir = dirname + '/Shell_Slices/'
files = os.listdir(radatadir)
nfiles = len(files)
files.sort()

print ('Reading ' + radatadir + files[-number_from_end] + ' ...')
a = Shell_Slices(radatadir + files[-number_from_end], '')

nshells = a.nr

#nrows = int(np.ceil(np.sqrt(nshells)))
#ncols = int(np.ceil(nshells/nrows))
#nplots = ncols*nrows  #some of these may be empty

# Make a list of the desired variable names
var_names = ['vr', 'vt', 'vp', 'vp_prime', 'vl', 'vz', \
                'omr', 'omt', 'omp', 'oml', 'omz',\
                'vsq', 'vhsq', 'vrsq',  'omsq', 's',\
                'vrs', 'vrvp', 'vrvt', 'vtvp',\
                'vlvp', 'vzvp', 'vlvz']

# For the variables that are directly output by Rayleigh
var_names_v = ['vr', 'vt', 'vp']
var_names_om = ['omr', 'omt', 'omp']
var_names_indexed = var_names_v + var_names_om
var_names_indexed.pop() # Nick's vort_phi is messed up!

if not oldra:
    var_indices = {\
        'vr'    :       1, 
        'vt'    :       2,
        'vp'    :       3,
        'omr'   :       301,
        'omt'   :       302,
        'omp'   :       303, 
        's'     :       501 }
else:
    var_indices = {\
        'vr'    :       1, 
        'vt'    :       2,
        'vp'    :       3,
        'omr'   :       103,
        'omt'   :       104,
        'omp'   :       105,
        's'     :       64 }


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

# Set up the lat/lon grid for the Basemap projections
xpixels = 1024  
ypixels = 1024

#Generate 1-D grids of latitude and longitude
nphi = a.nphi
nt = a.ntheta
cost = a.costheta
sint = a.sintheta
dlon = 360.0/(nphi)
dlat = 180.0/nt
lons = np.zeros(nphi)
for i in range(nphi):
    lons[i] = dlon*i-180.0

lats = np.zeros(nt)
for i in range(nt):
    lats[i] = 90.0-np.arccos(cost[i])*180.0/np.pi
    lats[i] = i*dlat-90

# Make a 2-D lat/lon grid (non-cyclic) for the orthographic projection
llons, llats = np.meshgrid(lons, lats)

# Create the plots. This first portion only handles the projection 
# and colorbars. Labeling comes further down.

for varname in ['vr']:
#for varname in ['vr', 'vsq']:
#for varname in ['vr']: # for debugging purposes
    posdef = False
    if (varname in ['vrsq', 'vhsq', 'vsq', 'omsq']):
        posdef = True
    # Make the directory to save the plots if it doesn't exist
    savedir = plotdir + 'sslice/' + varname + '/'
    if (not os.path.isdir(savedir)):
        os.makedirs(savedir)
#    for i_time in range(a.niter): # range over each time in the iteration
    for i_time in range(1): # for debugging purposes
        iter_loc = a.iters[i_time]

        for i_depth in range(nshells): # range over each depth
#        for i_depth in range(1): # for debugging purposes
            common_name = '_' + str(iter_loc).zfill(8) + '_depth' + str(i_depth).zfill(2) + '.png'
            savename_moll = 'moll' + common_name
            savename_ortho = 'ortho' + common_name
            print(savedir)
            print(savename_moll)
            print(savename_ortho)
            r_loc = a.radius[i_depth]
            ir = a.inds[i_depth]
            i_plot = i_depth + 1
            
            # Compute the location of the tangent cylinder at this depth
#            tangent_theta = np.arcsin(ri/r_loc)
#            tangent_theta_deg = 180./np.pi*tangent_theta
#            tangent_lat = 90. - tangent_theta_deg

            iter_loc = a.iters[i_time]

            vr_slice = a.vals[:, :, i_depth,\
                    a.lut[var_indices['vr']], i_time]
            vt_slice = a.vals[:, :, i_depth,\
                    a.lut[var_indices['vt']], i_time]
            vp_slice = a.vals[:, :, i_depth,\
                    a.lut[var_indices['vp']], i_time]
#            vp_slice_prime = vp_slice - (vp_av[:, ir]).reshape((1,nt))
            vp_slice_prime = vp_slice -\
                    (np.mean(vp_slice, axis=1)).reshape((nphi,1))
            s_slice = a.vals[:, :, i_depth,\
                    a.lut[var_indices['s']], i_time]

#            try:
#                omr_slice = a.vals[:, :, i_depth,\
#                        a.lut[var_indices['omr']], i_time]
#                omt_slice = a.vals[:, :, i_depth,\
#                        a.lut[var_indices['omt']], i_time]
#                omp_slice = a.vals[:, :, i_depth,\
#                        a.lut[var_indices['omp']], i_time]
#            except:
#                pass

#            try: # only plot if possible, then move on
#            if (True):
                # First, do the Mollweide projection
            print ('plotting ' + varname + ', depth ' + \
                str(i_depth).zfill(2) + ', iter ' + \
                str(iter_loc).zfill(8) + ', Mollweide  ...')
           
            if (varname == 'vr'):
                sslice = vr_slice/100. # measure velocities in m/s
            elif (varname == 'vt'):
                sslice = vt_slice/100.
            elif (varname == 'vp'):
                sslice = vp_slice/100.
            elif (varname == 'vp_prime'):
                sslice = vp_slice_prime/100.

             # for the cylindrical coordinate vectors, apply the 
             # appropriate rotation transformations
            elif (varname == 'vl'):
                sslice = vr_slice*sint_nd + vt_slice*cost_nd
                sslice /= 100.
            elif (varname == 'vz'):
                sslice = vr_slice*cost_nd - vt_slice*sint_nd
                sslice /= 100.
            
             # vorticities
            elif (varname == 'omr'):
                sslice = omr_slice
            elif (varname == 'omt'):
                sslice = omt_slice

             # cylindrical vorticities
            elif (varname == 'oml'):
                sslice = omr_slice*sint_nd + omt_slice*cost_nd
            elif (varname == 'omz'):
                sslice = omr_slice*cost_nd - omt_slice*sint_nd

            elif (varname == 'omp'):
                sslice = omp_slice - vt_slice + vt_slice/r_loc

             # velocities/vorticities-squared (vsq --> (m/s)^2)
            elif (varname == 'vrsq'):
                sslice = vr_slice**2
                sslice /= 1.0e4
            elif (varname == 'vhsq'):
                sslice = vt_slice**2 + vp_slice**2
                sslice /= 1.0e4
            elif (varname == 'vsq'):
                sslice = vr_slice**2 + vt_slice**2 + vp_slice**2
                sslice /= 1.0e4
            elif (varname == 'omsq'):
                omp = omp_slice - vt_slice + vt_slice/r_loc
                sslice = omp**2 + omr_slice**2 + omt_slice**2

             # the entropy fluctuation
            elif (varname == 's'):
#                    sslice = s_slice - s_spherical_mean[ir]
                sslice = s_slice - np.mean(s_slice)
            
            # velocity and entropy (self-) correlations
            elif (varname == 'vrs'):
                sslice = (s_slice - s_spherical_mean[ir])*vr_slice/100.
            elif (varname == 'vrvp'):
                sslice = vr_slice*vp_slice_prime/1.0e4
            elif (varname == 'vrvt'):
                sslice = vr_slice*vt_slice/1.0e4
            elif (varname == 'vtvp'):
                sslice = vt_slice*vp_slice_prime/1.0e4
            elif (varname == 'vlvp'):
                sslice = (vr_slice*sint_nd + vt_slice*cost_nd)*vp_slice_prime/1.0e4
            elif (varname == 'vzvp'): 
                sslice = (vr_slice*cost_nd - vt_slice*sint_nd)*vp_slice_prime/1.0e4
            elif (varname == 'vlvz'):
                sslice = (vr_slice*sint_nd + vt_slice*cost_nd)*\
                (vr_slice*cost_nd - vt_slice*sint_nd)/1.0e4
                
            sslice = np.transpose(sslice)

            # Make the individual plots. Start with Mollweide

            # Set up the cyclic grid for Mollweide projection
            sslice_cycl, lons_cycl = addcyclic(sslice, lons)

            # Convert to 2-D grids
            llons_cycl, llats_cycl = np.meshgrid(lons_cycl, lats)
            
            # Make the Mollweide plot
            fig = plt.figure(figsize=(6,3), dpi=300)
            m = Basemap(projection='moll', lon_0=0, resolution='l')

            # Get x-y projection points on the plane
            x, y = m(llons_cycl, llats_cycl)

            #View the data
            my_cmap = plt.cm.RdYlBu_r
            if (posdef):
                my_cmap = 'Greys'

            sig = np.std(sslice)
            my_max = 3*sig
            my_min = -3*sig

            # Saturation levels for a positive-definite quantity (like vsq)
            logslice = np.log(sslice)
            medlog = np.median(logslice)
            shiftlog = logslice - medlog
            minexp = medlog - 7*np.std(shiftlog[np.where(shiftlog < 0)].flatten())
            maxexp = medlog + 7*np.std(shiftlog[np.where(shiftlog > 0)].flatten())
            if not posdef:
                m.pcolormesh(x, y, sslice_cycl, cmap=my_cmap, vmin=my_min, vmax=my_max)
            else: 
                m.pcolormesh(x,y,sslice_cycl, cmap='Greens', norm=colors.LogNorm(vmin=np.exp(minexp),\
                        vmax=np.exp(maxexp)))

            # draw parallels and meridians. Draw two parallels at the tangent cylinder
#            m.drawparallels((-tangent_lat, tangent_lat))
            m.drawmeridians(np.arange(0.,420.,60.))
            m.drawmapboundary(fill_color='aqua')
            
            if not posdef:
                cbar = m.colorbar(format=ticker.FuncFormatter(fmt))
            else:
                cbar = m.colorbar()

            cbar.set_label(units[varname], rotation=0, labelpad = 15)
            cbar.ax.tick_params(labelsize=8)

#            percent_depth = (ro - r_loc)/d*100
#            plt.title(labels[varname] + ', depth = %.1f%%' %percent_depth)
            plt.tight_layout()
            plt.savefig(savedir + savename_moll)
            plt.close()

            # Then, do the orthographic projection
            print ('plotting ' + varname + ', depth ' + \
                str(i_depth).zfill(2) + ', iter ' + \
                str(iter_loc).zfill(8) + ', orthographic  ...')

            f1 = plt.figure(figsize=(5.3,4), dpi=300)
            ro = np.max(a.radius)
            m = Basemap(projection='ortho',lon_0=0,lat_0=20,resolution=None,rsphere=ro)

            topodat,x,y =\
            m.transform_scalar(sslice, lons, lats, xpixels, ypixels, returnxy=True, masked=True, order=1)

            shrink_distance = ro - r_loc
            shrink_factor = ro/r_loc
            x = x/shrink_factor + shrink_distance
            y = y/shrink_factor + shrink_distance 

            #View the data
            if not posdef:
                m.pcolormesh(x, y, topodat, cmap=my_cmap, vmin=my_min, vmax=my_max)
            else: 
                m.pcolormesh(x,y, topodat, cmap='Greens', norm=colors.LogNorm(vmin=np.exp(minexp),\
                        vmax=np.exp(maxexp)))

            # add a color bar
            if not posdef:
                cbar = m.colorbar(format=ticker.FuncFormatter(fmt))
            else:
                cbar = m.colorbar()
            cbar.set_label(units[varname], rotation=0, labelpad = 15)
            cbar.ax.tick_params(labelsize=8)

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
#                for tlat in [-tangent_lat, tangent_lat]:
#                    linex, liney = m(lons, np.ones_like(lons)*tlat)
#                    # must cut off all points on the far side of the sphere
#                    linex = linex[np.where(np.abs(linex) < 1.e20)]
#                    liney = liney[np.where(np.abs(liney) < 1.e20)]
#
#                    linex = linex/shrink_factor + shrink_distance
#                    liney = liney/shrink_factor + shrink_distance

#                    m.plot(linex, liney,'k',linewidth=.75)

#                plt.title(labels[varname] + ', depth = %.1f%%' %percent_depth)
            plt.tight_layout()
            plt.subplots_adjust(left=0)
            plt.savefig(savedir + savename_ortho)
            plt.close()
#            except:
#                pass
