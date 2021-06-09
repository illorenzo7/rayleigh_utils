# Author: Loren Matilsky
# Date created: 03/02/2019
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
from common import *
from plotcommon import axis_range
from tl_util import plot_tlr
from cla_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# get data
if 'the_file' in clas: 
    the_file = clas['the_file']
else:
    the_file = get_widest_range_file(clas0['datadir'], 'timelat_b')

# Read in the time-latitude data (dictionary form)
print ('Getting time-latitude trace from ' + the_file)
di = get_dict(the_file)
vals = di['vals']
times = di['times']
iters = di['iters']
rvals_avail = di['samplevals']
qvals_avail = np.array(di['qvals'])

# baseline time unit
iter1, iter2 = get_iters_from_file(the_file)
time_unit, time_label, rotation, simple_label = get_time_unit(dirname)

br_index = np.argmin(np.abs(qvals - 801))
bt_index = np.argmin(np.abs(qvals - 802))
bp_index = np.argmin(np.abs(qvals - 803))

# determine desired levels to plot
if not 'irvals' in clas:
    if not 'rvals' in clas:
        irvals = np.array([0]) # just plot the top radius by default
    else: # get irvals from rvals
        rvals = clas['rvals']
        if rvals == 'all':
            irvals = np.arange(a.nr)
        else:
            irvals = np.zeros_like(rvals, dtype='int')
            for i in range(len(rvals)):
                irvals[i] = np.argmin(np.abs(rvals_avail - rvals[i]))
else:
    irvals = clas['irvals']

# Get raw traces of br, btheta, bphi
br = vals[:, :, :, br_index]
bt = vals[:, :, :, bt_index]
bp = vals[:, :, :, bp_index]

# Normalize the time 
times /= time_unit

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

nrow = 3
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

# field units and labels
units = r'$\rm{G}$'
labels = [r'$\langle B_r\rangle$', r'$\langle B_\theta\rangle$',\
        r'$\langle B_\phi\rangle$']

# get grid info
di_grid = get_grid_info(dirname)

# Loop over the desired radii and save plots
for i in range(len(irvals)):
    ir = irvals[i]
    rval = rvals_avail[ir]
    print('plotting r/rsun = %0.3f (ir = %02i)' %(rval, ir))
    br_loc = br[:, :, ir]
    bt_loc = bt[:, :, ir]
    bp_loc = bp[:, :, ir]
    
    # Make appropriate file name to save
    savename = 'timelat_B_' + ('%08i_%08i_' %(iter1, iter2)) +\
        ('rval%0.3f' %rval) + '.png'

    if 'minmax' in clas:
        minmax = clas['minmax']
    else:
        minmax = None
    if minmax is None:
        minmax_br = None
        minmax_bt = None
        minmax_bp = None
    else:
        if len(minmax) == 2:
            minmax_br = minmax
            minmax_bt = minmax
            minmax_bp = minmax
        elif len(minmax) == 6:
            minmax_br = minmax[0], minmax[1]
            minmax_bt = minmax[2], minmax[3]
            minmax_bp = minmax[4], minmax[5]
        else:
            print ("error: minmax must have length 2 or 6, as in")
            print ("-minmax '-5e3 5e3 -1e4 1e4 -5e4 5e4'")
            print ("exiting")
            sys.exit()

    # Create figure with  3 panels in a row (time-radius plots of
    #       br, btheta, and bphi)
    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax1 = fig.add_axes((margin_left, 1. - margin_top -\
            subplot_height - 0*(subplot_height + margin_y),\
            subplot_width, subplot_height))
    ax2 = fig.add_axes((margin_left, 1. - margin_top -\
            subplot_height - 1*(subplot_height + margin_y),\
            subplot_width, subplot_height))
    ax3 = fig.add_axes((margin_left, 1. - margin_top -\
            subplot_height - 2*(subplot_height + margin_y),\
            subplot_width, subplot_height))

    # Plot evolution of each (zonally averaged) field component
    plot_tlr(br_loc, times, di_grid['tt_lat'], fig, ax1,\
            minmax=minmax_br)
    plot_tlr(bt_loc, times, di_grid['tt_lat'], fig, ax2,\
            minmax=minmax_bt)
    plot_tlr(bp_loc, times, di_grid['tt_lat'], fig, ax3,\
            minmax=minmax_bp)

    # Label with the field components
    for irow in range(nrow):
        label = labels[irow]
        fig.text(margin_left + 0.5*margin_x, 1. - margin_top - \
                0.5*margin_y - irow*(subplot_height + margin_y), label,\
                va='top', ha='left', fontsize=14,\
                bbox=dict(facecolor='white'))

    # Turn the x tick labels off for the top two strips
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])

    # Label x (time) axis
    ax3.set_xlabel('time (' + time_label + ')')
    # Label y-axis (latitude in degrees)
    ax2.set_ylabel('latitude (deg)')

    # Put some useful information on the title
    title = dirname_stripped + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %rval)
    if 'navg' in clas:
        navg = clas['navg']
        averaging_time = (times[-1] - times[0])/len(times)*navg
        title += '     ' + ('t_avg = %.1f Prot' %averaging_time)
    else: 
        title += '     t_avg = none'
    ax1.set_title(title)

    # Save the plot
    if clas0['saveplot']:
        # save the figure
        plotdir = my_mkdir(clas0['plotdir'] + 'timelat/')
        print ('Saving the time-latitude plot at ')
        print (plotdir + savename)
        print ("=======================================")
        #plt.savefig(plotdir + savename, dpi=200)

    # Show the plot if only plotting at one latitude
    #if clas0['showplot'] and len(irvals) == 1:
    plt.show()
    plt.close()
