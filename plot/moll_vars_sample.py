import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['co'])
sys.path.append(os.environ['rapp'])
from common import get_file_lists, strip_dirname, rsun
from sslice_util import plot_moll, get_sslice, axis_range
from rayleigh_diagnostics import Shell_Slices
from get_parameter import get_parameter
from translate_times import translate_times
from varprops import texlabels

# Get command line arguments
dirname = sys.argv[1]
radatadir = dirname + '/Shell_Slices/'

file_list, int_file_list, nfiles = get_file_lists(radatadir)

iiter = nfiles - 1 # by default plot the last iteration
ir = 0 # by default plot just below the surface
rval = None # can also find ir by finding the closest point
clon = 0.

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-ir'):
        ir = int(args[i+1])
    elif arg == '-rval':
        rval = float(args[i+1])
    elif (arg == '-clat'):
        clat = float(args[i+1])
    elif (arg == '-clon'):
        clon = float(args[i+1])
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

iter_val = int_file_list[iiter]
fname = file_list[iiter]

# Read in desired shell slice
a = Shell_Slices(radatadir + fname, '')

# Find desired radius (by default ir=0--near outer surface)
if not rval is None:
    ir = np.argmin(np.abs(a.radius/rsun - rval))
rval = a.radius[ir] # in any case, this is the actual rvalue we get

# Create the plot using subplot axes
# Offset axes slightly (at the end) to deal with annoying white space cutoff
fig_width_inches = 12.

# General parameters for main axis/color bar
margin_bottom_inches = 1.
margin_top_inches = 1.
margin_inches = 1/8

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

# Loop over depths and make plots
plotdir = dirname + '/plots/moll/vars_sample/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Loop over vars and make plots
varlist = ['vr_prime', 'vt_prime', 'vp_prime', 's_prime',\
        'p_prime', 's_prime_sph']
magnetism = get_parameter(dirname, 'magnetism')
if magnetism:
    varlist.extend(['br', 'bt', 'bp'])

for varname in varlist: 
    savename = 'moll_iter' + fname + ('_rval%0.3f_' %(rval/rsun)) +\
            varname  + '.png'
    print('Plotting moll: ' + varname + (', r/rsun = %0.3f (ir = %02i), '\
            %(rval/rsun, ir)) + 'iter ' + fname + ' ...')
    vals = get_sslice(a, varname, dirname=dirname)
    field = vals[:, :, ir]

    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax = fig.add_axes([margin_x, margin_bottom, subplot_width,\
            subplot_height])

    plot_moll(field, a.costheta, fig=fig, ax=ax, varname=varname,\
            clon=clon) 

    # Make title
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x    
    
    varlabel = texlabels[varname]
    title = varlabel + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %(rval/rsun)) +\
            '     ' + ('iter = ' + fname)    
    fig.text(ax_center_x, ax_ymax + 0.02*ax_delta_y, title,\
         verticalalignment='bottom', horizontalalignment='center',\
         fontsize=14, **csfont)   
    
    fig.text(margin_x + 0.5*subplot_width, 1. - 0.5*margin_top,\
        strip_dirname(dirname), ha='center', va='bottom', **csfont,\
        fontsize=14)

    plt.savefig(plotdir + savename, dpi=300)
    plt.close()
