# Created 03/05/2020
# Shows the meridional slice associated with a variable without saving
# the plot
# Similar to moll_show
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from common import get_file_lists, strip_dirname, rsun
from plotcommon import axis_range
from translate_times import translate_times
from azav_util import plot_azav
from rayleigh_diagnostics import Meridional_Slices, AZ_Avgs, Shell_Avgs
from get_merslice import get_merslice
from get_parameter import get_parameter
from varprops import texlabels, texunits
from time_scales import compute_Prot, compute_tdt

# Get command line arguments
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Data with Shell_Slices
radatadir = dirname + '/Meridional_Slices/'
file_list, int_file_list, nfiles = get_file_lists(radatadir)

minmax = None
symlog = False
logscale = False
iiter = nfiles - 1 # by default plot the last iteration
iphi = -1 # by default plot the meridian closest to 0 longitude
varname = 'vr' # by default plot the radial velocity
plotcontours = True
plotlatlines = True
use_az = True
use_sh = True

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-iphi':
        iphi = int(args[i+1])
    elif arg == '-var':
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
vals = get_merslice(mer, varname, dirname=dirname, sh=sh, az=az)

# Get local time (in seconds)
t_loc = mer.time[0]

field = vals[iphi, :, :]
lonval = mer.phi[iphi]*180./np.pi

# Figure dimensions
subplot_width_inches = 2.5
subplot_height_inches = 5.
margin_inches = 1./8.
margin_top_inches = 1. # larger top margin to make room for titles

fig_width_inches = subplot_width_inches + 2*margin_inches
fig_height_inches = subplot_height_inches + margin_top_inches +\
        margin_inches

margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

# create axes
units = texunits[varname]
texlabel = texlabels[varname]

fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
ax = fig.add_axes((margin_x, margin_y, subplot_width, subplot_height))
plot_azav (field, mer.radius, mer.costheta, fig=fig, ax=ax, units=units,\
        minmax=minmax, plotlatlines=plotlatlines, plotcontours=plotcontours)

# Make title + label diff. rot. contrast and no. contours
# Label what time it is
if rotation:
    time_string = ('t = %.1f ' %(t_loc/time_unit)) + time_label
else:
    time_string = ('t = %.3f ' %(t_loc/time_unit)) + time_label

fsize = 12.
line_spacing_inches = 1./4.
space = line_spacing_inches/fig_height_inches
fig.text(margin_x, 1 - space, dirname_stripped, ha='left',\
        va='bottom', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 2*space,\
         time_string, ha='left', va='bottom', fontsize=fsize,\
         **csfont)
fig.text(margin_x, 1 - 3*space, texlabel +\
        (r'$\ \ \ \ \ \phi = %03.1f^\circ$' %lonval),\
         ha='left', va='bottom', fontsize=fsize, **csfont)
plt.show()