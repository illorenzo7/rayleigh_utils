# Author: Loren Matilsky
# Date created: 08/18/2019
#import matplotlib as mpl
#mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['raco'])
from common import *
from plotcommon import axis_range

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]

# Data and plot directories
datadir = dirname + '/data/'
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Get global rotation rate; this script fails for non-rotating models
angular_velocity = get_parameter(dirname, 'angular_velocity')
Prot = 2*np.pi/angular_velocity
Om0_nhz = angular_velocity/(2*np.pi)*1e9
print('Om0: ', Om0_nhz, ' nHz')

# more defaults
minmax = None # this may be true at first
minmax_wasnone = True # This always will be true unless user specifies
            # values through -minmax
# By default (if tminmax is None) will plot over whole time trace
tminmax = None
saveplot = None # turned off by default if saving one figure, can change
# with -save option
showplot = False # will only show if plotting one figure
labelbytime = False # by default label by first/last iteration number

rvals = 'all'  # by default, plot all available time-lat levels
    # user specifies another choice via -rvals '[val1] [val2] ... [valn]'
    # where 'vals' have dimensional units in cm: 4.8e10, 5e10, etc.
irvals = None # user can also specify -irvals '2 3 9', etc.

# By default, plot longitudinal B field
qval = 803

Om_subtract = None # by default, do not subtract any differential rotation
# if nonzero, subtract a CONSTANT DR

# Default no. rotations to plot per inch (vertically)
rpi = 100.

# Get command-line arguments
plotdir = None

args = sys.argv[2:]
nargs = len(args)
tag = ''
plottimes = None
clat = 10. 
dlat = 0.
the_file = None
for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
        minmax_wasnone = False
    elif arg == '-usefile':
        the_file = args[i+1]
        the_file = the_file.split('/')[-1]
    elif arg == '-rvals':
        strings = args[i+1].split()
        rvals = []
        for j in range(len(strings)):
            rvals.append(float(strings[j]))
    elif arg == '-irvals':
        irvals = []
        strings = args[i+1].split()
        for j in range(len(strings)):
            irvals.append(int(strings[j]))
    elif arg == '-tminmax':
        tminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-save':
        saveplot = True
    elif arg == '-tlabel':
        labelbytime = True
    elif arg == '-rpi': # rpi = rotations per inch
        rpi = float(args[i+1])
    elif arg == '-qval':
        qval = int(args[i+1])
    elif arg == '-tag':
        tag = args[i+1] + '_'
    elif arg == '-clat':
        clat = float(args[i+1])
    elif arg == '-dlat':
        dlat = float(args[i+1])
    elif arg == '-om': # enter this value in nHz, in the LAB frame
        Om_subtract = float(args[i+1]) - Om0_nhz
        Om_subtract *= (2*np.pi/1e9) # convert nHz --> rad s^-1
    elif arg == '-times':
        strings = args[i+1].split()
        plottimes = []
        for string in strings:
            plottimes.append(float(string))

# Get plot directory and create if not already there
if plotdir is None:
    plotdir = dirname + '/plots/time-lon' + tag + '/'
    if labelbytime:
        plotdir = dirname + '/plots/time-lon' + tag + '_tlabel/'
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)

# Find the time-longitude file(s) the data directory, for clat and dlat. 
# If there are 
# multiple, by default choose the one with widest range in the trace.
if the_file is None:
    if clat >= 0.:
        hemisphere = 'N'
    else:
        hemisphere = 'S'

    the_file = get_widest_range_file(datadir, 'time-longitude_clat' + hemisphere +\
            '%02.0f_dlat%03.0f' %(np.abs(clat), dlat))

# Read in the time-longitude data (dictionary form)
print ('Getting time-longitude trace from ' + datadir + the_file)
di = get_dict(datadir + the_file)

vals = di['vals']
times = di['times']
iters = di['iters']
lut = di['lut']
qvals = di['qvals']
rr = di['rr']
ri = di['ri']; ro = di['ro']; shell_depth = ro - ri
clat = di['clat']
dlat = di['dlat']
irvals_avail = di['rinds']
rvals_avail = rr[irvals_avail]
lons = di['lons']
nphi = di['nphi']

niter = di['niter']
nr = di['nr']
nrvals = di['nrvals']
nq = di['nq']

iter1 = di['iter1']
iter2 = di['iter2']

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

# determine desired levels to plot
if irvals is None:
    if rvals == 'all':
        irvals = np.arange(len(irvals_avail))
    else:
        irvals = []
        for rval in rvals:
            ir = np.argmin(np.abs(rvals_avail - rval))
            irvals.append(ir)
if saveplot is None:
    if len(irvals) == 1:
        saveplot = False
    else:
        saveplot = True
if len(irvals) == 1:
    showplot = True

# Get raw trace of "qval"
if tminmax is None:
    it1, it2 = 0, niter - 1
else:
    it1 = np.argmin(np.abs(times/Prot - tminmax[0]))
    it2 = np.argmin(np.abs(times/Prot - tminmax[1]))

q_index = np.argmin(np.abs(qvals - qval))

quant = vals[it1:it2+1, :, :, q_index]
times = times[it1:it2+1]/Prot
t1, t2 = times[0], times[-1] # These begin times and end times
        # will be used for labeling the plots

# Make meshgrid of time/radius
lons_2d, times_2d = np.meshgrid(lons, times, indexing='ij')

# Subtract DR, if desired
if not Om_subtract is None:
    phi_deflections = (times*Prot*Om_subtract*180./np.pi) % 360.
    for it in range(len(times)):
        phi_deflection = phi_deflections[it]
        nroll = int(phi_deflection/360*nphi)
        quant[it, :, :] = np.roll(quant[it, :, :], -nroll, axis=0)

# Loop over the desired radii and save plots
for i in range(len(irvals)):
    ir = irvals[i]
    rval = rvals_avail[ir]/rsun 
    print('plotting r/rsun = %0.3f (ir = %02i)' %(rval, ir))
    
    quant_loc = quant[:, :, ir]

    # Make appropriate file name to save
    if labelbytime:
        savename = dirname_stripped + '_time-lon_' +\
                ('Prot%05.0f-to-%05.0f_clat%s%02.0f_dlat%02.0f_'\
                %(t1, t2, hemisphere, np.abs(clat), dlat)) +\
                ('qv%04i_' %qval) + ('rval%0.3f' %rval) + '.png'
    else:
        savename = dirname_stripped + '_time-lon_' +\
                ('%08i_%08i_clat%s%02.0f_dlat%02.0f_'\
                %(iter1, iter2, hemisphere, np.abs(clat), dlat)) +\
                ('qv%04i_' %qval) + ('rval%0.3f' %rval) + '.png'

    if minmax is None or minmax_wasnone:
        std_quant = np.std(quant_loc)
        minmax = -3.*std_quant, 3.*std_quant
        minmax_wasnone = True
     
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
    ax.set_xlabel(r'$\rm{longitude}\ (^\circ)$', fontsize=12)

    # Label y (time) axis
    timeunit = r'$$'
    ylabel = r'$\rm{time}\ (P_{\rm{rot}})}$'
#    ylabel = r'$\rm{time}\ (\longleftarrow)}$'
    ax.set_ylabel(ylabel, **csfont, fontsize=12)

    ax.set_xlim((lons[0], lons[-1]))

    # Label the plots by B_r, B_theta, B_phi
    ax_xmin, ax_xmax = ax.get_xlim()
    ax_Dx = ax_xmax - ax_xmin
    ax.text(ax_xmin + 1.01*ax_Dx, 0.,  'qval = %i' %qval, **csfont)

    # Put some useful information on the title
    title = dirname_stripped + (', iq = %i' %qval) + '\n' +\
            (r'$r/R_\odot\ =\ %0.3f$' %rval) + '\n' +\
            (r'$\lambda_c=%02.0f^\circ\ \ \ \ \ \Delta\lambda=%02.0f^\circ$'\
            %(clat, dlat)) 
    ax.set_title(title, **csfont)

    # Get ticks everywhere
    plt.sca(ax)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    # Put grid of white lines on plot
    ax.grid(color='k', linestyle='-', linewidth=1., alpha=0.3)
    ax.grid(which='minor', color='k', linestyle='-', linewidth=0.5, alpha=0.3)

    #Save the plot
    if saveplot:
        print ('Saving the time-latitude plot at ' + plotdir +\
                savename)
        plt.savefig(plotdir + savename, dpi=300)

    # Show the plot if only plotting at one latitude
    if showplot:
        plt.show()
    plt.close()
