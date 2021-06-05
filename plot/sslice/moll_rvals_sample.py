import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
sys.path.append(os.environ['rapp'])
from common import *
from sslice_util import plot_moll
from get_sslice import get_sslice
from rayleigh_diagnostics import Shell_Slices
from get_slice import get_slice, get_label

# Get CLAs
args = sys.argv 
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Data with Shell_Slices
radatadir = dirname + '/Shell_Slices/'
args = sys.argv[2:]
nargs = len(args)
file_list, int_file_list, nfiles = get_file_lists(radatadir, args)

minmax = None
iiter = nfiles - 1 # by default plot the last iteration
varname = 'vr' # by default plot the radial velocity
clon = 0.
rvals = []

plotdir = None

for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-var' or arg == '-qval':
        varname = args[i+1]
    elif arg == '-iter':
        desired_iter = int(args[i+1])
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-sec':
        time = float(args[i+1])
        di_trans = translate_times(time, dirname, translate_from='sec')
        desired_iter = di_trans['val_iter']
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-day':
        time = float(args[i+1])
        di_trans = translate_times(time, dirname, translate_from='day')
        desired_iter = di_trans['val_iter']
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-prot':
        time = float(args[i+1])
        di_trans = translate_times(time, dirname, translate_from='prot')
        desired_iter = di_trans['val_iter']
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-clon':
        clon = float(args[i+1])
    elif arg == '-depths':
        strings = args[i+1].split()
        for st in strings:
            rval = ro - float(st)*d
            rvals.append(rval)
    elif arg == '-depthscz':
        rm = domain_bounds[1]
        dcz = ro - rm
        strings = args[i+1].split()
        for st in strings:
            rval = ro - float(st)*dcz
            rvals.append(rval)
    elif arg == '-depthsrz':
        rm = domain_bounds[1]
        drz = rm - ri
        strings = args[i+1].split()
        for st in strings:
            rval = rm - float(st)*drz
            rvals.append(rval)
    elif arg == '-rvals':
        rvals = []
        strings = args[i+1].split()
        for st in strings:
            rval = float(st)*rsun
            rvals.append(rval)
    elif arg == '-rvalscm':
        rvals = []
        strings = args[i+1].split()
        for st in strings:
            rval = float(st)
            rvals.append(rval)
    elif arg == '-rrange':
        r1 = float(args[i+1])
        r2 = float(args[i+2])
        n = int(args[i+3])
        rvals = np.linspace(r1, r2, n)*rsun

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

if plotdir is None:
    plotdir = dirname + '/plots/moll/rvals_sample/'
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)

# File name to read
iter_val = int_file_list[iiter]
fname = file_list[iiter]

# Read in desired shell slice
a = Shell_Slices(radatadir + fname, '')

# Get local time (in seconds)
t_loc = a.time[0]

# Create the plot using subplot axes
# Offset axes slightly (at the end) to deal with annoying white space cutoff
fig_width_inches = 6.

# General parameters for main axis/color bar
margin_bottom_inches = 1./2.
margin_top_inches = 5./8.
margin_inches = 1./8.

subplot_width_inches = fig_width_inches - 2*margin_inches
subplot_height_inches = 0.5*subplot_width_inches
fig_height_inches = margin_bottom_inches + subplot_height_inches +\
    margin_top_inches
fig_aspect = fig_height_inches/fig_width_inches

# "Non-dimensional" figure parameters
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches

subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

# Loop over rvals and make plots
if rvals == []:
    irvals = np.arange(a.nr) # loop over everything
else:
    irvals = np.zeros_like(rvals, dtype='int')
    for i in range(len(rvals)):
        irvals[i] = np.argmin(np.abs(a.radius - rvals[i]))

for ir in irvals:
    rval = a.radius[ir]/rsun
    vals = get_sslice(a, varname, dirname=dirname)
    field = vals[:, :, ir]

    savename = 'moll_' + varname + '_iter' + fname +\
            ('_rval%0.3f' %rval) + '.png'
    # Display at terminal what we are plotting
    print('Plotting moll: ' + varname + (', r/rsun = %0.3f (ir = %02i), '\
            %(rval, ir)) + 'iter ' + fname)

    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax = fig.add_axes([margin_x, margin_bottom, subplot_width,\
            subplot_height])
    
    plot_moll(field, a.costheta, fig=fig, ax=ax,\
            minmax=minmax, clon=clon) 

    # Make title
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x    
    
    if rotation:
        time_string = ('t = %.1f ' %(t_loc/time_unit)) + time_label +\
                ' (1 ' + time_label + (' = %.2f days)'\
                %(time_unit/86400.))
    else:
        time_string = ('t = %.3f ' %(t_loc/time_unit)) + time_label +\
                ' (1 ' + time_label + (' = %.1f days)'\
                %(time_unit/86400.))
    varlabel = texlabels.get(varname, 'qval = ' + varname)

    title = dirname_stripped +\
        '\n' + r'$\rm{Mollweide}$' + '     '  + time_string +\
        '\n' + varlabel + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %rval)
    fig.text(ax_center_x, ax_ymax + 0.02*ax_delta_y, title,\
         verticalalignment='bottom', horizontalalignment='center',\
         fontsize=10, **csfont)   
    plt.savefig(plotdir + savename, dpi=300)
    plt.close()
