# Author: Loren Matilsky
# Created: 06/11/2018
####################################################################################################
#
#  Shell-Slice (Shell_Slices) plotting example
#  - Reads in a single Shell_Slice file.
#  - Plots orthographic projections of the shell slice at each depth available
#  - By default plots vr, though other variables may be indicated with the
#  - "-var" option (.e.g., -var s, for entropy)
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
sys.path.append('./rayleigh/plot')
from common import get_widest_range_file, strip_dirname
from sslice_util import plot_ortho, labels, units, cbar_axes_from_area

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'

if not os.path.isdir(plotdir):
    os.makedirs(plotdir)
    
# Read command-line arguments (CLAs)
number_from_end = 1 # By default, plot the last slice
args = sys.argv[2:]
nargs = len(args)
varname = 'vr'
use_old_lut = False
use_specific_iter = False
for i in range(nargs):
    arg = args[i]
    if (arg == '-n'):
        number_from_end = int(args[i+1])
    elif (arg == '-plotdir'):
        plotdir = dirname + '/' + args[i+1] + '/'
    elif (arg == '-var'):
        varname = args[i+1]
    elif (arg == '-old'):
        use_old_lut = True
    elif (arg == '-iter'):
        use_specific_iter = True
        iter_to_use = int(args[i+1])

# Get grid info
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nt, nr = len(tt), len(rr)
nphi = 2*nt
cost_nd = cost.reshape((1, nt))
sint_nd = sint.reshape((1, nt))

# Get mean entropy
savg_file = get_widest_range_file(datadir, 's_spherical_mean')
s_spherical_mean = np.load(datadir + savg_file)

# Get mean velocities
vavg_file = get_widest_range_file(datadir, 'vavg')
vavg = np.load(datadir + vavg_file)
vr_av, vt_av, vp_av = vavg


# Make a list of the desired variable names
var_names = ['vr', 'vt', 'vp', 'vp_prime', 'vl', 'vz', \
                'omr', 'omt', 'omp', 'oml', 'omz',\
                'vsq', 'vhsq', 'vrsq',  'omsq', 's',\
                'vrs', 'vrvp', 'vrvt', 'vtvp',\
                'vlvp', 'vzvp', 'vlvz']

# Get a list of all the shellslices
radatadir = dirname + '/Shell_Slices/'
files = os.listdir(radatadir)
nfiles = len(files)
files.sort()

fname = files[-number_from_end]
if (use_specific_iter):
    fname = str(iter_to_use).zfill(8)

if (use_old_lut):
    from diagnostic_reading import ShellSlice
    reading_function = ShellSlice
else:
    from rayleigh_diagnostics import Shell_Slices
    reading_function = Shell_Slices
    
print ('Reading ' + radatadir + fname + ' ...')
a = reading_function(radatadir + files[-number_from_end], '')
iter_loc = a.iters[0]
nshells = a.nr

# initialize indices, depending on whether "-old" was indicated
qindex_vr = 1
qindex_vt = 2
qindex_vp = 3

qindex_omr = 301
qindex_omt = 302
qindex_omp = 303

qindex_s = 501

if (use_old_lut):
    qindex_omr = 103
    qindex_omt = 104
    qindex_omp = 105
    
    qindex_s = 64

vr_slice = a.vals[:, :, :, a.lut[qindex_vr], 0]
vt_slice = a.vals[:, :, :, a.lut[qindex_vt], 0]
vp_slice = a.vals[:, :, :, a.lut[qindex_vp], 0]
vp_slice_prime = vp_slice - (vp_av[:, a.inds]).reshape((1, nt, nshells))
s_slice = a.vals[:, :, :, a.lut[qindex_s], 0]

try: # I didn't always output the vorticities
    omr_slice = a.vals[:, :, :, a.lut[qindex_omr], 0]
    omt_slice = a.vals[:, :, :, a.lut[qindex_omt], 0]
    omp_slice = a.vals[:, :, :, a.lut[qindex_omp], 0]
except:
    pass

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
    sslice = omp_slice - vt_slice + vt_slice/(a.radius).reshape((1, 1, nshells))

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
    omp = omp_slice - vt_slice + vt_slice/(a.radius).reshape((1, 1, nshells))
    sslice = omp**2 + omr_slice**2 + omt_slice**2

 # the entropy fluctuation
elif (varname == 's'):
    sslice = s_slice - (s_spherical_mean[a.inds]).reshape((1, 1, nshells))

# velocity and entropy (self-) correlations
elif (varname == 'vrs'):
    sslice = (s_slice - (s_spherical_mean[a.inds]).reshape((1, 1, nshells)))*\
    vr_slice/100.
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

#Generate 1-D grids of latitude and longitude
dlon = 360.0/nphi
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

# Set up the actual figure from scratch
fig_width_inches = 3.5 # TOTAL figure width, in inches (i.e., 8x11.5 paper with 1/2-inch
                # margins and two text columns ~3.5 inches in width)
margin_inches = 1/8 # margin width in inches (for both x and y) and horizontally 
            # in between figures
margin_top_inches = 1/8

nrow = 1
ncol = 1

colorbar_area_inches = 1/2  # width of the area to which the whole colorbar 
    # (including labels) belongs
            
subplot_width_inches = (fig_width_inches - (ncol + 1)*margin_inches)/ncol

subplot_height_inches = subplot_width_inches # Each subplot should have an
    # aspect ratio of y/x = 1/1 for an orthographic projection.
fig_height_inches = margin_top_inches +\
    nrow*(subplot_height_inches + margin_inches + colorbar_area_inches)

fig_aspect = fig_height_inches/fig_width_inches

# "Margin" in "figure units"; figure units extend from 0 to 1 in BOTH 
# directions, so unitless dimensions of margin will be different in x and y
# to force an equal physical margin
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_left = margin_inches/fig_width_inches
margin_top = margin_top_inches/fig_height_inches
cbar_area = colorbar_area_inches/fig_height_inches

# Location of "main axis" (axis on figure accommodating the margins)
# in figure units
main_axis_left = margin_left
main_axis_bottom = cbar_area
main_axis_width = 1 - margin_left - margin_x
main_axis_height = 1 - margin_top - cbar_area

# Subplot dimensions in figure units
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches               

cbar_axes = cbar_axes_from_area(main_axis_left,\
    0, subplot_width, cbar_area, fig_aspect)
            
posdef = False
if (varname in ['vrsq', 'vhsq', 'vsq', 'omsq']):
    posdef = True
# Make the directory to save the plots if it doesn't exist
savedir = plotdir + 'sslice/all_depths/' + varname + '/'
if (not os.path.isdir(savedir)):
    os.makedirs(savedir)

for i_depth in range(nshells): # range over each depth
#for i_depth in range(1): # for debugging purposes
    sslice_depth = np.transpose(sslice[:, :, i_depth])
    print ('plotting ' + varname + ', depth ' + \
           str(i_depth).zfill(3) + ', iter ' + \
           str(iter_loc).zfill(8) + ', Orthographic projection ...')
    savename_ortho = dirname_stripped + '_ortho_' + varname + '_' +\
        str(iter_loc).zfill(8) + '_depth' + str(i_depth).zfill(3) + '.png'
    savename_moll = dirname_stripped + '_moll_' + varname + '_' +\
        str(iter_loc).zfill(8) + '_depth' + str(i_depth).zfill(3) + '.png'
    r_loc = a.radius[i_depth]
    ir = a.inds[i_depth]
    
    # Set up axes to add a subplot, if it's possible (sometimes have to deal 
        # with NaNs at the extreme points like depth000)
    try:
        fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
        
        ax = fig.add_axes((main_axis_left, main_axis_bottom,\
                            subplot_width, subplot_height))
        cax = fig.add_axes(cbar_axes, fig_aspect)
        plot_ortho(fig, ax, sslice_depth, lats, lons, a.radius[i_depth], ri, ro,\
                   cax=cax, cbar_units=units[varname], posdef=posdef)
        plt.title(labels[varname] + ', ' + r'$\ r/r_o=%.3f$' %(r_loc/ro))
        print ('Saving ' + savedir + savename_ortho + ' ...')
        plt.savefig(savedir + savename_ortho, dpi=300)
        plt.close()
    except:
        pass