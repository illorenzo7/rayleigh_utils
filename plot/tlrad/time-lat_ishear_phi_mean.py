# Author: Loren Matilsky
# Date created: 03/02/2019
import matplotlib as mpl
from matplotlib import ticker
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
from common import *
from plotcommon import axis_range
from tl_util import plot_tl

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]

# Data and plot directories
datadir = dirname + '/data/'
nosave = False
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Find the time/latitude file(s) the data directory. If there are 
# multiple, by default choose the one with widest range in the trace.
# We need both time-latitude regular (for the B field)
# and time-latitude induction (for the induction terms)
the_file = get_widest_range_file(datadir, 'time-latitude_shear_terms')

# more defaults
minmax = None
xminmax = None
xmin = None
xmax = None
saveplot = None # turned off by default if saving one figure, can change
# with -save option
showplot = False # only show if plotting one figure
labelbytime = False # by default label by first/last iteration number
# not first/last time
rvals = 'all'  # by default, plot all available time-lat levels
    # user specifies another choice via -rvals '[val1] [val2] ... [valn]'
    # where 'vals' have dimensional units in cm: 4.8e10, 5e10, etc.
irvals = None # user can also specify -irvals '2 3 9', etc.
navg = 1 # by default average over 1 AZ_Avgs instance (no average)
# for navg > 1, a "sliding average" will be used.
tag = '' # optional way to tag save directory
lats = [0.]
plottimes = None
phi_deriv = False

# Get command-line arguments
plotdir = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-minmax':
        strings = args[i+1].split()
        minmax = []
        for st in strings:
            minmax.append(float(st))
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
    elif arg == '-navg':
        navg = int(args[i+1])
        if navg % 2 == 0:
            print ("Please don't enter even values for navg!")
            print ("Replacing navg = %i with navg = %i" %(navg, navg + 1))
            navg += 1
    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-xmin':
        xmin = float(args[i+1])
    elif arg == '-xmax':
        xmax = float(args[i+1])
    elif arg == '-save':
        saveplot = True
    elif arg == '-tlabel':
        labelbytime = True
    elif arg == '-tag':
        tag = '_' + args[i+1]
    elif arg == '-lats':
        strings = args[i+1].split()
        for string in strings:
            lats.append(float(string))
    elif arg == '-times':
        strings = args[i+1].split()
        plottimes = []
        for string in strings:
            plottimes.append(float(string))
    elif arg == '-dp':
        phi_deriv = True

# Get plot directory and create if not already there
plotdir = dirname + '/plots/time-lat' + tag + '_ishear_phi_mean/'
if labelbytime:
    plotdir = dirname + '/plots/time-lat' + tag + '_ishear_phi_mean_tlabel/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Read in the time-latitude data (dictionary form)
# In this part, we really must assume that the user saved
# the separate time-latitude (induction) files wisely, so they all
# have the same values for "times" and "lats", depths, etc.
print ('Getting time-latitude induction trace from ' + datadir +\
       the_file)
di = get_dict(datadir + the_file)

vals = di['vals']
vals = di['vals']

times = di['times']
iters = di['iters']
rr = di['rr']
irvals_avail = di['rinds']
rvals_avail = rr[irvals_avail]
tt_lat = di['tt_lat']
rinds = di['rinds'] # radial locations sampled for the trace
ntheta = di['ntheta']
rvals_sampled = rr[rinds]/rsun

niter = di['niter']
nr = di['nr']
ntimes = di['niter']

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

if plotdir is None:
    plotdir = dirname + '/plots/'
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)
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

# shear terms: total, B_r*dv/dr, B_t*dv/dt, B_p*dv/dp, 2 curvature terms
npp = 6
if phi_deriv:
    npp += 1
ind_off = 2*(5 + npp)
terms = []
for i in range(6):
    terms.append(vals[:, :, :, ind_off + i])
if phi_deriv:
    terms.insert(4, vals[:, :, :, ind_off + 6])

# field units and labels
units = r'$\rm{G\ s^{-1}}$'
labels = [r'$[\left\langle \mathbf{B}\right\rangle\cdot\nabla\left\langle\mathbf{v}\right\rangle]_\phi$',\
r'$\left\langle B_r\right\rangle\left\langle\frac{\partial v_\phi}{\partial r}\right\rangle$',\
r'$\frac{1}{r}\left\langle B_\theta\right\rangle\left\langle\frac{\partial v_\phi}{\partial\theta}\right\rangle$',\
r'$\frac{1}{r\sin\theta}\left\langle B_\phi\right\rangle\left\langle\frac{\partial v_\phi}{\partial\phi}\right\rangle$',\
r'$\frac{1}{r}\langle B_\phi\rangle \langle v_r\rangle$',\
r'$\frac{\cot\theta}{r}\langle B_\phi\rangle \langle v_\theta\rangle}$']
if phi_deriv:
    labels.insert(4, r'$\frac{1}{r\sin\theta}\left\langle B_\phi\right\rangle\left\langle\frac{\partial v_\phi}{\partial\phi}\right\rangle$' + ' exact')

# Normalize the time 
times /= time_unit

# Make meshgrid of time/latitude
# Take into account if user specified xmin, xmax
if xminmax is None:
    xminmax = np.min(times), np.max(times)
# Change JUST xmin or xmax, if desired
if not xmin is None:
    xminmax = xmin, xminmax[1]
if not xmax is None:
    xminmax = xminmax[0], xmax

it1 = np.argmin(np.abs(times - xminmax[0]))
it2 = np.argmin(np.abs(times - xminmax[1]))
t1, t2 = times[it1], times[it2] # These begin times and end times
        # will be used for labeling the plots

# set figure dimensions
subplot_width_inches = 6.5
margin_inches = 1./4.
margin_bottom_inches = 1./2. # space for x-axis label
margin_top_inches = 1./2.
margin_left_inches = 5./8. # space for latitude label
margin_right_inches = 0.9

fig_width_inches = subplot_width_inches + margin_right_inches +\
        margin_left_inches
subplot_height_inches = 2.0

nrow = len(terms)
fig_height_inches = nrow*subplot_height_inches +\
        (nrow - 1)*margin_inches + margin_bottom_inches +\
        margin_top_inches

margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches
margin_left = margin_left_inches/fig_width_inches
margin_bottom = margin_bottom_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches

# Loop over the desired radii and save plots
for i in range(len(irvals)):
    ir = irvals[i]
    rval = rvals_avail[ir]/rsun 
    print('plotting r/rsun = %0.3f (ir = %02i)' %(rval, ir))
    terms_loc = []
    for j in range(nrow):
        terms_loc.append(terms[j][:, :, ir])
   
    # Make appropriate file name to save
    if labelbytime:
        savename = dirname_stripped + '_time-lat_ishear_phi_mean_' +\
                ('Prot%05.0f-to-%05.0f_' %(t1, t2)) +\
            ('rval%0.3f' %rval) + '.png'
    else:
        savename = dirname_stripped + '_time-lat_ishear_phi_mean_' +\
                ('%08i_%08i_' %(iter1, iter2)) +\
            ('rval%0.3f' %rval) + '.png'

    mins_and_maxes = []
    if minmax is None:
        for j in range(nrow):
            mins_and_maxes.append(None)
    else:
        if len(minmax) == 2:
            for j in range(nrow):
                mins_and_maxes.append(minmax)
        elif len(minmax) == 2*nrow:
            for j in range(nrow):
                mins_and_maxes.append((minmax[2*j], minmax[2*j + 1]))
        else:
            print ("error: minmax must have length 2 or %i, as in"%(2*nrow))
            print ("-minmax '-5e-2 5e-2'")
            print ("your minmax has length %i" %(len(minmax)))
            print ("exiting")
            sys.exit()

    # make plots subplot by subplot
    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    axs = []
    for j in range(nrow):
        # Create figure with  3 panels in a row (time-radius plots of
        #       br, btheta, and btheta)
        axs.append(fig.add_axes((margin_left, 1. - margin_top -\
                subplot_height - j*(subplot_height + margin_y),\
                subplot_width, subplot_height)))

        # Plot evolution of each (zonally averaged) field component
        plot_tl(terms_loc[j], times, tt_lat, fig=fig, ax=axs[j], navg=navg,\
                minmax=mins_and_maxes[j], units=units, xminmax=xminmax,\
                yvals=lats, plottimes=plottimes)

        # Label each subplot
        fig.text(margin_left + 0.5*margin_x, 1. - margin_top - \
                0.5*margin_y - j*(subplot_height + margin_y), labels[j],\
                va='top', ha='left', fontsize=14,\
                bbox=dict(facecolor='white'))

    # Turn the x tick labels off for the top strips
    for ax in axs[:-1]:
        ax.set_xticklabels([])

    # Label x (time) axis
    axs[-1].set_xlabel('time (' + time_label + ')', **csfont)
    # Label y-axis (latitude in degrees)
    axs[nrow//2].set_ylabel('latitude (deg)', **csfont)

    # Put some useful information on the title
    averaging_time = (times[-1] - times[0])/niter*navg
    title = dirname_stripped + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %rval)
    if navg > 1:
        title += '     ' + ('t_avg = %.1f Prot' %averaging_time)
    else:
        title += '     t_avg = none'
    axs[0].set_title(title, **csfont)

    # Save the plot
    if saveplot:
        print ('Saving the time-latitude plot at ')
        print (plotdir + savename)
        print ("=======================================")
        plt.savefig(plotdir + savename, dpi=200)

    # Show the plot if only plotting at one latitude
    if showplot:
        plt.show()
    plt.close()