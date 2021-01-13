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
from sslice_util import plot_ortho
from get_sslice import get_sslice
from rayleigh_diagnostics import Shell_Slices
from varprops import texlabels
from get_eq import get_eq

# Get command line arguments
dirname = sys.argv[1]
radatadir = dirname + '/Shell_Slices/'

file_list, int_file_list, nfiles = get_file_lists(radatadir)

minmax = None
iiter = nfiles - 1 # by default plot the last iteration
varname = 'vr' # by default plot the radial velocity

clon = 0
clat = 20

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-clat':
        clat = float(args[i+1])
    elif arg == '-clon':
        clon = float(args[i+1])        
    elif arg == '-var':
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
        
iter_val = int_file_list[iiter]
fname = file_list[iiter]

# Read in desired shell slice
a = Shell_Slices(radatadir + fname, '')

# We also need the radius, which we can get from the reference state
eq = get_eq(dirname)
radius = eq.radius

# Create the plot template
fig_width_inches = 6.

# General parameters for main axis/color bar
margin_bottom_inches = 1.
margin_top_inches = 1.
margin_inches = 1/8

subplot_width_inches = fig_width_inches - 2*margin_inches
subplot_height_inches = subplot_width_inches
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

# Make plot (sub-)directory if it doesn't already exist
plotdir = dirname + '/plots/ortho/rvals_sample/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Loop over rvals and make plots   
for ir in range(a.nr):
    rval = a.radius[ir]
    vals = get_sslice(a, varname, dirname=dirname)
    field = vals[:, :, ir]
    
    savename = 'ortho_iter' + fname     
    if clat >= 0:
        hemisphere = 'N'
    else:
        hemisphere = 'S'
    savename += ('_clat' + hemisphere + '%02.0f_' %(np.abs(clat))) 
    savename += (varname + ('_rval%0.3f.png' %(rval/rsun)))

    print('Plotting ortho: ' + varname + (', r/rsun = %0.3f (ir = %02i), '\
            %(rval/rsun, ir)) + 'iter ' + fname +\
            ', clat = %02.0f' %clat + ' ...')
    
    # Make axes and plot the orthographic projection
    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax = fig.add_axes([margin_x, margin_bottom, subplot_width,\
            subplot_height])
    plot_ortho(field, radius, a.costheta, fig=fig, ax=ax, ir=a.inds[ir],\
        clon=clon, clat=clat, varname=varname)  
    
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
