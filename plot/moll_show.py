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
from sslice_util import plot_moll
from rayleigh_diagnostics import Shell_Slices
from get_sslice import get_sslice
from varprops import texlabels

# Get command line arguments
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Data with Shell_Slices
radatadir = dirname + '/Shell_Slices/'

file_list, int_file_list, nfiles = get_file_lists(radatadir)

minmax = None
symlog = False
logscale = False
iiter = nfiles - 1 # by default plot the last iteration
ir = 0 # by default plot just below the surface
rval = None # can also find ir by finding the closest point
            # to a local radius divided by rsun
varlist = ['vr'] # by default plot the radial velocity
clon = 0
saveplot = False
ncol = 2

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-ir':
        ir = int(args[i+1])
    elif arg == '-rval':
        rval = float(args[i+1])
    elif arg == '-var' or arg == '-qval':
        varlist = args[i+1].split()
    elif arg == '-symlog':
        symlog = True
    elif arg == '-log':
        logscale = True
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
    elif arg == '-save':
        saveplot = True
        tag = args[i+1] + '_'
    elif arg == '-ncol':
        ncol = int(args[i+1])

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

iter_val = int_file_list[iiter]
fname = file_list[iiter]

# Read in desired shell slice
a = Shell_Slices(radatadir + fname, '')

# figure dimensions
nplots = len(varlist)
if ncol > nplots:
    ncol = nplots # (no point in having empty columns)
nrow = np.int(np.ceil(nplots/ncol))
fig_width_inches = 6.*ncol

# General parameters for main axis/color bar
margin_bottom_inches = 1./2.
margin_top_inches = 3./4. # margin for big title
margin_subplot_top_inches = 1/4 # margin to accommodate just subplot titles
margin_inches = 1./8.

subplot_width_inches = (fig_width_inches - (ncol + 1)*margin_inches)\
    /ncol
subplot_height_inches = 0.5*subplot_width_inches
fig_height_inches = margin_top_inches + nrow*(margin_bottom_inches +\
    subplot_height_inches + margin_subplot_top_inches)
fig_aspect = fig_height_inches/fig_width_inches

# "Non-dimensional" figure parameters
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_subplot_top = margin_subplot_top_inches/fig_height_inches

subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

# Main figure
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

# loop over desired vars and make plots
for iplot in range(nplots):
    varname = varlist[iplot]
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1. - margin_top - margin_subplot_top -\
            subplot_height - (iplot//ncol)*(subplot_height +\
            margin_subplot_top + margin_bottom)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))

    vals = get_sslice(a, varname, dirname=dirname)

    # Get local time (in seconds)
    t_loc = a.time[0]

    # Find desired radius (by default ir=0--near outer surface)
    if not rval is None:
        ir = np.argmin(np.abs(a.radius/rsun - rval))
    field = vals[:, :, ir]
    rval = a.radius[ir]/rsun # in any case, this is the actual rvalue we get

    # Display at terminal what we are plotting
    print('Plotting moll: ' + varname + (', r/rsun = %0.3f (ir = %02i), '\
            %(rval, ir)) + 'iter ' + fname)

    plot_moll(field, a.costheta, fig=fig, ax=ax, minmax=minmax, clon=clon,\
            varname=varname, symlog=symlog, logscale=logscale) 

    # label the subplot
    varlabel = texlabels.get(varname, 'qval = ' + varname)
    ax.set_title(varlabel, verticalalignment='bottom', **csfont)

    # Make title
    if iplot == 0:
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

        title = dirname_stripped +\
            '\n' + r'$\rm{Mollweide}$' + '     '  + time_string +\
            '\n' + (r'$r/R_\odot\ =\ %0.3f$' %rval)
        fig.text(ax_center_x, ax_ymax + 0.02*ax_delta_y +\
                margin_subplot_top, title,\
             verticalalignment='bottom', horizontalalignment='center',\
             fontsize=10, **csfont)   
if saveplot:
    plotdir = dirname + '/plots/moll/'
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)
    savefile = plotdir + 'moll_' + tag + fname + '.png'
    print ('Saving plot at ' + savefile)
    plt.savefig(savefile, dpi=300)
# always show
plt.show()   
