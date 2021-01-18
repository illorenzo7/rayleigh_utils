# Created 04/23/2020
# Plots and saves meridional slices associated for a range of variables at 
# a given longitude and time
# Similar to moll_vars_sample
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
iphi = -1 # by default plot the meridian closest to 0 longitude
lonval = None # can also find iphi by finding the closest point
            # to a local longitude 
iiter = nfiles - 1 # by default plot the last iteration
plotcontours = True
plotlatlines = True
use_az = True
use_sh = True
varlist = None
rbcz = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-minmaxrz':
        minmaxrz = float(args[i+1]), float(args[i+2])
    elif arg == '-rbcz':
        rbcz = float(args[i+1])
    elif arg == '-var':
        varname = args[i+1]
    elif arg == '-symlog':
        symlog = True
    elif arg == '-nocontour':
        plotcontours = False
    elif arg == '-nolat':
        plotlatlines = True
    elif arg == '-noaz':
        use_az = False
    elif arg == '-nosh':
        use_sh = False
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
    elif arg == '-vars':
        varlist = args[i+1].split()

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

# Check if there are corresponding Shell_Avgs or AZ_Avgs files if needed
if use_az:
    try:
        az = AZ_Avgs(dirname + '/AZ_Avgs/' + fname, '')
        print ("read AZ_Avgs/" + fname)
    except:
        print("No file AZ_Avgs/" + fname)
        print("setting az = None")
        az = None
else:
    az = None

if use_sh:
    try:
        sh = Shell_Avgs(dirname + '/Shell_Avgs/' + fname, '')
        print ("read Shell_Avgs/" + fname)
    except:
        print("No file Shell_Avgs/" + fname)
        print("setting sh = None")
        sh = None
else:
    sh = None

# Read in desired meridional slice
print ("read Meridional_Slices/" + fname)
mer = Meridional_Slices(radatadir + fname, '')

# Get local time (in seconds)
t_loc = mer.time[0]

# Get the longitude
# Find desired longitude value (by default iphi=-1--near 0 deg)
if not lonval is None:
    iphi = np.argmin(np.abs(mer.phi*180./np.pi - lonval))
lonval = mer.phi[iphi]*180./np.pi
# in any case, this is the actual lonvalue we get

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

# Make plot (sub-)directory if it doesn't already exist
plotdir = dirname + '/plots/mer/vars_sample/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Loop over vars and make plots
if varlist is None:
    varlist = ['vr_prime', 'vt_prime', 'vp_prime', 'omr_prime',\
            'omt_prime', 'omp_prime', 's_prime', 'p_prime', 's_prime_sph']
    magnetism = get_parameter(dirname, 'magnetism')
    if magnetism:
        varlist.extend(['br', 'bt', 'bp'])

for varname in varlist: 
    # Get Tex units and variable label
    units = texunits.get(varname, 'cgs')
    texlabel = texlabels.get(varname, varname)

    vals = get_merslice(mer, varname, dirname=dirname, sh=sh, az=az)
    field = vals[iphi, :, :]
    savename = 'mer_iter' + fname + ('_lonval%05.1f_' %lonval) +\
            varname  + '.png'

    # Display at terminal what we are plotting
    print('Plotting mer: ' + varname + (', lon = %.1f (ir = %02i), '\
            %(lonval, iphi)) + 'iter ' + fname)

    # Create axes
    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax = fig.add_axes((margin_x, margin_bottom, subplot_width,\
            subplot_height))
    plot_azav (field, mer.radius, mer.costheta, fig=fig, ax=ax,\
            units=units, minmax=minmax, minmaxrz=minmaxrz,\
            plotlatlines=plotlatlines, plotcontours=plotcontours, rbcz=rbcz)

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
    fig.text(margin_x, 1. - margin_y, dirname_stripped_title, ha='left',\
            va='top', fontsize=fsize, **csfont)
    fig.text(margin_x, 1. - margin_y - 2*line_height, 'Meridional Slice',\
            ha='left', va='top', fontsize=fsize, **csfont)
    fig.text(margin_x, 1. - margin_y - 3*line_height,\
             time_string, ha='left', va='top', fontsize=fsize,\
             **csfont)
    fig.text(margin_x, 1 - margin_y - 5*line_height, texlabel +\
            (r'$\ \ \ \ \ \phi = %03.1f^\circ$' %lonval),\
             ha='left', va='top', fontsize=fsize, **csfont)
    plt.savefig(plotdir + savename, dpi=300)
    plt.close()
