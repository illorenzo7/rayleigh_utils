# Author: Loren Matilsky
# Date created: 08/18/2019
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['co'])
from common import get_file_lists, get_widest_range_file, strip_dirname,\
        rsun, get_dict
from plotcommon import axis_range
from get_parameter import get_parameter

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]

# Data and plot directories
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)
dirname_stripped = strip_dirname(dirname)

# Find the time-longitude file(s) the data directory. If there are 
# multiple, by default choose the one with widest range in the trace.
time_longitude_file = get_widest_range_file(datadir, 'time-longitude')
print('the file', time_longitude_file)

# more defaults
minmax = None
xminmax = None
saveplot = True
showplot = True # will only show if plotting one figure

desired_rvals = [0.83] # by default, plot time-radius diagram for fields 
    # mid-CZ (units of solar radius)

# By default, plot longitudinal B field
qval = 803

# By default (if tminmax is None) will plot over whole time trace
tminmax = None

# Default no. rotations to plot per inch (vertically)
rpi = 100.
# Get command-line arguments
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-usefile':
        time_longitude_file = args[i+1]
        time_longitude_file = time_latitude_file.split('/')[-1]
    elif arg == '-rvals':
        string_desired_rvals = args[i+1].split()
        if string_desired_rvals == ['all']:
            desired_rvals = 'all'
        else:
            desired_rvals = []
            for j in range(len(string_desired_rvals)):
                desired_rvals.append(float(string_desired_rvals[j]))
    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-nosave':
        saveplot = False
    elif arg == '-noshow':
        showplot = False
    elif arg == '-rpi': # rpi = rotations per inch
        rpi = float(args[i+1])
    elif arg == '-qval':
        qval = int(args[i+1])

# Read in the time-longitude data (dictionary form)
print ('Getting time-longitude trace from ' + datadir +\
       time_longitude_file + ' ...')
di = get_dict(datadir + time_longitude_file)

vals = di['vals']
times = di['times']
iters = di['iters']
lut = di['lut']
rr = di['rr']
ri = di['ri']; ro = di['ro']; shell_depth = ro - ri
clat = di['clat']
dlat = di['dlat']
rinds = di['rinds'] # radial locations sampled for the trace
rvals_sampled = rr[rinds]/rsun
lons = di['lons']
nphi = di['nphi']

niter = di['niter']
nr = di['nr']
nrvals = di['ndepths']
nq = di['nq']

iter1 = di['iter1']
iter2 = di['iter2']

# Get global rotation rate; this script fails for non-rotating models
angular_velocity = get_parameter(dirname, 'angular_velocity')
Prot = 2*np.pi/angular_velocity

# Get radial indices of values we want to plot
i_desiredrvals = []
rvals_to_plot = []
if desired_rvals == 'all':
    i_desiredrvals = np.arange(len(rinds))
    rvals_to_plot = rvals_sampled
else:
    i_desiredrvals = []
    rvals_to_plot = []
    for desired_rval in desired_rvals:
        i_desiredrval = np.argmin(np.abs(rvals_sampled - desired_rval))
        i_desiredrvals.append(i_desiredrval)
        rvals_to_plot.append(rvals_sampled[i_desiredrval])

# Get raw trace of "qval"
if tminmax is None:
    it1, it2 = 0, niter - 1
else:
    it1 = np.argmin(np.abs(times/Prot - tminmax[0]))
    it2 = np.argmin(np.abs(times/Prot - tminmax[1]))

quant = vals[it1:it2+1, :, :, lut[qval]]
times = times[it1:it2+1]/Prot
t1, t2 = times[0], times[-1] # These begin times and end times
        # will be used for labeling the plots

# Make meshgrid of time/radius
lons_2d, times_2d = np.meshgrid(lons, times, indexing='ij')

# Loop over the desired radii and save plots
for i in range(len(i_desiredrvals)):
    i_desiredrval = i_desiredrvals[i]
    rval_to_plot = rvals_to_plot[i]
    
    quant_loc = quant[:, :, i_desiredrval]
    
    # Make appropriate file name to save
    savename = dirname_stripped + '_time-longitude_%i_' +\
            ('Prot%05.0f-to-%05.0f_clat%02.0f_dlat%02.0f_' %(t1, t2, clat, dlat)) +\
        ('rval%0.3f' %rval_to_plot) + '.png'

    if minmax is None:
        std_quant = np.std(quant_loc)
        minmax = -3.*std_quant, 3.*std_quant
     
    # Create figure as a vertical strip, taking into account the desired
    # rpi
    fig_width_inches = 5.
    fig_height_inches = (times[-1] - times[0])/rpi
    print('figsize = ', fig_width_inches, fig_height_inches)
    fig, ax = plt.subplots(figsize=(fig_width_inches, fig_height_inches)) 

    im = ax.pcolormesh(lons_2d, times_2d, quant_loc.T,\
            vmin=minmax[0], vmax=minmax[1], cmap='RdYlBu_r')
    plt.gca().invert_yaxis()

    # Put colorbar next to plot 
    # First make room and then find location of subplot
    plt.subplots_adjust(left=0.2, right=0.8, wspace=0.03, top=0.9)

    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x

    cbar_width = 0.15*(1 - ax_xmax)
    cbar_height = 0.4*ax_delta_y

    cbar_left = ax_xmax + 0.1*(1 - ax_xmax)
    cbar_bottom = ax_ymin + 0.5*ax_delta_y - 0.5*cbar_height
    cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))
    cax.set_title('cgs', **csfont)
    plt.colorbar(im, cax=cax)

    # Set x label
    ax.set_xlabel(r'$\rm{longitude}\ (^\circ)$')

    # Label y (time) axis
    timeunit = r'$$'
    ylabel = r'$\rm{time}\ (P_{\rm{rot}})}$'
#    ylabel = r'$\rm{time}\ (\longleftarrow)}$'
    ax.set_ylabel(ylabel, **csfont)

    ax.set_xlim((lons[0], lons[-1]))

    # Label the plots by B_r, B_theta, B_phi
    ax_xmin, ax_xmax = ax.get_xlim()
    ax_Dx = ax_xmax - ax_xmin
    ax.text(ax_xmin + 1.01*ax_Dx, 0.,  'qval = %i' %qval, **csfont)

    # Put some useful information on the title
    title = dirname_stripped + (', iq = %i' %qval) + '\n' +\
            (r'$r/R_\odot\ =\ %0.3f$' %rval_to_plot) + '\n' +\
            (r'$\lambda_c=%02.0f^\circ\ \ \ \ \ \Delta\lambda=%02.0f^\circ$' %(clat, dlat)) 
    ax.set_title(title, **csfont)

    # Get ticks everywhere
    plt.sca(ax)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    # Put grid of white lines on plot
    ax.grid(color='k', linestyle='-', linewidth=0.5, alpha=0.5)

    # Save the plot
    if saveplot:
        print ('Saving the time-latitude plot at ' + plotdir +\
                savename + ' ...')
#        plt.savefig(plotdir + savename, dpi=200)

    # Show the plot if only plotting at one latitude
    if len(rvals_to_plot) == 1 and showplot:
        plt.show()
    plt.close()
