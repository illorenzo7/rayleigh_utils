# Created 04/23/2020
# Plots and saves meridional slice associated with a given variable at 
# a given longitude in a time series
# Similar to moll_times_sample
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from common import *
from plotcommon import axis_range
from azav_util import plot_azav
from rayleigh_diagnostics import Meridional_Slices, AZ_Avgs, Shell_Avgs
from get_merslice import get_merslice
from varprops import texlabels, texunits

# Get command line arguments
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Split dirname_stripped into two lines if it is very long
if len(dirname_stripped) > 25:
    dirname_stripped_title = dirname_stripped[:25] + '\n' +\
            dirname_stripped[25:]
else:
    dirname_stripped_title = dirname_stripped

# Data with Meridional_Slices
radatadir = dirname + '/Meridional_Slices/'
file_list, int_file_list, nfiles = get_file_lists(radatadir)

minmax = None
minmaxrz = None
symlog = False
logscale = False
iphi = -1 # by default plot the meridian closest to 0 longitude
lonval = None # can also find iphi by finding the closest point
            # to a local longitude 
varname = 'vr' # by default plot the radial velocity
plotcontours = True
plotlatlines = True
use_az = True
use_sh = True
rbcz = None

plotdir = None

args = sys.argv[2:]
nargs = len(args)

the_tuple = get_desired_range(int_file_list, args)
if the_tuple is None: # By default plot the last 10 Shell_Slices
    index_first, index_last = nfiles - 11, nfiles - 1  
else:
    index_first, index_last = the_tuple

for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-minmaxrz':
        minmaxrz = float(args[i+1]), float(args[i+2])
    elif arg == '-rbcz':
        rbcz = float(args[i+1])
    elif arg == '-var' or arg == '-qval':
        varname = args[i+1]
    elif arg == '-symlog':
        symlog = True
    elif arg == '-log':
        logscale = True
    elif arg == '-nocontour':
        plotcontours = False
    elif arg == '-nolat':
        plotlatlines = True
    elif arg == '-noaz':
        use_az = False
    elif arg == '-nosh':
        use_sh = False
    elif arg == '-iphi':
        iphi = int(args[i+1])

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

# Figure dimensions
subplot_width_inches = 2.5
subplot_height_inches = 5.
margin_inches = 1./8.
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)
margin_top_inches = 1.75 # larger top margin to make room for titles

fig_width_inches = subplot_width_inches + 2*margin_inches
fig_height_inches = subplot_height_inches + margin_top_inches +\
        margin_bottom_inches

margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

# Get Tex units and variable label
units = texunits.get(varname, 'cgs')
texlabel = texlabels.get(varname, 'qval = ' + varname)

# Make plot (sub-)directory if it doesn't already exist
plotdir = dirname + '/plots/mer/times_sample/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Loop over each desired iteration and make plots
for i in range(index_first, index_last + 1):
    # Read in desired meridional slice
    fname = file_list[i]
    print ("read Meridional_Slices/" + fname)
    mer = Meridional_Slices(radatadir + fname, '')
    # Check if there are corresponding Shell_Avgs or AZ_Avgs files if needed
    if use_az and varname[-5:] == 'prime':
        try:
            az = AZ_Avgs(dirname + '/AZ_Avgs/' + fname, '')
            print ("read AZ_Avgs/" + fname)
        except:
            print("No file AZ_Avgs/" + fname)
            print("setting az = None")
            az = None
    else:
        az = None

    if use_sh and 'prime_sph' in varname:
        try:
            sh = Shell_Avgs(dirname + '/Shell_Avgs/' + fname, '')
            print ("read Shell_Avgs/" + fname)
        except:
            print("No file Shell_Avgs/" + fname)
            print("setting sh = None")
            sh = None
    else:
        sh = None

    # Find desired longitude value (by default iphi=-1--near 0 deg)
    if not lonval is None:
        iphi = np.argmin(np.abs(mer.phi*180./np.pi - lonval))
    lonval = mer.phi[iphi]*180./np.pi
    # in any case, this is the actual lonvalue we get

    # Loop over the slices in the file 
    for j in range(mer.niter):
        vals = get_merslice(mer, varname, dirname=dirname, sh=sh, az=az,\
                j=j)
        # Get local time (in seconds)
        t_loc = mer.time[j]
        iter_loc = mer.iters[j]

        # Get field associated with particular time/longitude
        field = vals[iphi, :, :]
        savename = 'mer_' + varname + ('_lonval%05.1f' %lonval) + '_iter' +\
                    str(iter_loc).zfill(8) + '.png'

        # Display at terminal what we are plotting
        print('Plotting mer: ' + varname + (', lon = %.1f (iphi = %02i), '\
                %(lonval, iphi)) + 'iter ' + str(iter_loc).zfill(8))

        # Create axes
        fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
        ax = fig.add_axes((margin_x, margin_bottom, subplot_width,\
                subplot_height))
        plot_azav (field, mer.radius, mer.costheta, fig=fig, ax=ax,\
                units=units, minmax=minmax, minmaxrz=minmaxrz,\
                plotlatlines=plotlatlines, plotcontours=plotcontours,\
                rbcz=rbcz)

        # Make title
        if rotation:
            time_string = ('t = %.1f ' %(t_loc/time_unit)) + time_label +\
                    '\n (1 ' + time_label + (' = %.2f days)'\
                    %(time_unit/86400.))
        else:
            time_string = ('t = %.3f ' %(t_loc/time_unit)) + time_label +\
                    '\n (1 ' + time_label + (' = %.1f days)'\
                    %(time_unit/86400.))

        fsize = 12.
        line_height = 1./4./fig_height_inches
        fig.text(margin_x, 1. - margin_y, dirname_stripped_title,\
                ha='left', va='top', fontsize=fsize, **csfont)
        fig.text(margin_x, 1. - margin_y - 2*line_height,\
                'Meridional Slice', ha='left', va='top',\
                fontsize=fsize, **csfont)
        fig.text(margin_x, 1. - margin_y - 3*line_height,\
                 time_string, ha='left', va='top', fontsize=fsize,\
                 **csfont)
        fig.text(margin_x, 1 - margin_y - 5*line_height, texlabel +\
                (r'$\ \ \ \ \ \phi = %03.1f^\circ$' %lonval),\
                 ha='left', va='top', fontsize=fsize, **csfont)
        plt.savefig(plotdir + savename, dpi=300)
        plt.close()
