# Author: Loren Matilsky
# Created: 04/06/2019
# This script plots mean magnetic field components in the meridional plane 
# for the Rayleigh run directory indicated by [dirname]. 
# Uses an instantaneous AZ_Avgs file unless told otherwise
# Saves plot in
# [dirname]_B_azav_[iter].npy
# or [dirname]_B_azav_[first iter]_[last iter].npy if a time average was 
# specified

import numpy as np
import pickle
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
from binormalized_cbar import MidpointNormalize
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from azav_util import plot_azav
from common import get_widest_range_file, strip_dirname, get_file_lists,\
        get_desired_range
from get_parameter import get_parameter
from time_scales import compute_Prot, compute_tdt
from translate_times import translate_times
from rayleigh_diagnostics import AZ_Avgs

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

radatadir = dirname + '/AZ_Avgs/'

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Set defaults
mins = None
maxes = None
save = True
plotcontours = True

# Read in CLAs (if any) to change default variable ranges and other options
args = sys.argv[2:]
nargs = len(args)

if '-range' in args or '-n' in args or '-centerrange' in args or\
        '-all' in args or '-iter' in args:
    index_first, index_last = get_desired_range(int_file_list, args)
else:
    index_first, index_last = nfiles - 1, nfiles - 1  
    # By default, don't average over any files;
    # just plot the last file

# Get the time range in sec
t1 = translate_times(int_file_list[index_first], dirname,\
        translate_from='iter')['val_sec']
t2 = translate_times(int_file_list[index_last], dirname,\
        translate_from='iter')['val_sec']

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

# Change other defaults
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        mins = float(args[i+1]), float(args[i+3]), float(args[i+5])
        maxes = float(args[i+2]), float(args[i+4]), float(args[i+6])
    elif arg == '-nosave':
        saveplot = False
    elif arg == '-nocontour':
        plotcontours = False

# Initialize empty "vals" array for the time average
az0 = AZ_Avgs(radatadir + file_list[index_first], '')
rr = az0.radius
ri, ro = np.min(rr), np.max(rr)
d = ro - ri
rr_depth = (ro - rr)/d
rr_height = (rr - ri)/d
sint = az0.sintheta
cost = az0.costheta
tt = np.arccos(cost)
tt_lat = (np.pi/2 - tt)*180/np.pi
nr = az0.nr
nt = az0.ntheta

# compute some derivative quantities for the grid
tt_2d, rr_2d = np.meshgrid(tt, rr, indexing='ij')
sint_2d = np.sin(tt_2d); cost_2d = np.cos(tt_2d)
xx = rr_2d*sint_2d
zz = rr_2d*cost_2d

count = 0
iter1, iter2 = int_file_list[index_first], int_file_list[index_last]

vals = np.zeros_like(az0.vals[:, :, :, 0])

# Average over the relevant data range, summing everything and then dividing
#   by the number of "slices" added at the end
print ('Considering AZ_Avgs files %s through %s for the average ...'\
       %(file_list[index_first], file_list[index_last]))
for i in range(index_first, index_last + 1):
    print ('Adding AZ_Avgs/%s to the average ...' %file_list[i])
    if i == index_first:
        az = az0
    else:   
        az = AZ_Avgs(radatadir + file_list[i], '')

    local_ntimes = az.niter
    times_to_loop = np.arange(local_ntimes)
    if index_last == index_first:
        times_to_loop = np.array([-1])
    for j in times_to_loop:
        vals += az.vals[:, :, :, j]
        count += 1

vals /= count
print ('Averaged over %i AZ_Avgs slice(s) ...' %count)

br = vals[:, :, az0.lut[801]]
bt = vals[:, :, az0.lut[802]]
bp = vals[:, :, az0.lut[803]]

if mins is None and maxes is None:
    nstd = 5.
    mins = -nstd*np.std(br), -nstd*np.std(bt), -nstd*np.std(bp)
    maxes = nstd*np.std(br), nstd*np.std(bt), nstd*np.std(bp)

# Set up the actual figure from scratch
fig_width_inches = 7. # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1./8. # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_top_inches = 1.5 # wider top margin to accommodate subplot titles AND metadata
margin_subplot_top_inches = 1. # margin to accommodate just subplot titles
ncol = 3 # put three plots per row
nrow = 1

subplot_width_inches = (fig_width_inches - (ncol + 1)*margin_inches)/ncol
    # Make the subplot width so that ncol subplots fit together side-by-side
    # with margins in between them and at the left and right.
subplot_height_inches = 2*subplot_width_inches # Each subplot should have an
    # aspect ratio of y/x = 2/1 to accommodate meridional planes. 
fig_height_inches = nrow*subplot_height_inches + margin_top_inches +\
    (nrow - 1)*margin_subplot_top_inches + margin_inches 
    # Room for titles on each row and a regular margin on the bottom
fig_aspect = fig_height_inches/fig_width_inches

# "Margin" in "figure units"; figure units extend from 0 to 1 in BOTH 
# directions, so unitless dimensions of margin will be different in x and y
# to force an equal physical margin
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_subplot_top = margin_subplot_top_inches/fig_height_inches

# Subplot dimensions in figure units
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

field_components = [br, bt, bp]
titles = [r'$B_r$', r'$B_\theta$', r'$B_\phi$']
units = r'$\rm{G}$'

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

for iplot in range(3):
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1 - margin_top - subplot_height - \
            (iplot//ncol)*(subplot_height + margin_subplot_top)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    plot_azav (field_components[iplot], rr, cost, fig=fig, ax=ax,\
           units=units, minmax = (mins[iplot], maxes[iplot]),\
           plotcontours=plotcontours)
    ax.set_title(titles[iplot], verticalalignment='bottom', **csfont)

# Label averaging interval
if rotation:
    time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
            + time_label + '\n' + (r'$\ (\Delta t = %.1f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'
else:
    time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.3f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'

# Put some metadata in upper left
fsize = 12
line_height = 1./4./fig_height_inches
fig.text(margin_x, 1 - margin_y, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - margin_y - line_height,\
        'Magnetic field (zonally averaged)', ha='left', va='top',\
        fontsize=fsize, **csfont)
fig.text(margin_x, 1 - margin_y - 2*line_height, time_string,\
         ha='left', va='top', fontsize=fsize, **csfont)

savename = dirname_stripped + '_B_azav_' + str(iter1).zfill(8) + '_' +\
        str(iter2).zfill(8) + '.png'
if save:
    plt.savefig(plotdir + savename, dpi=100)
    print ("Saving azimuthal average of B field in " + plotdir + savename + ' ...')
plt.show()
