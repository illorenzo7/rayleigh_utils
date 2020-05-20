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
from azav_util import plot_azav, streamfunction
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
showplot = True
minmax = None
saveplot = True
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

# Change other defaults
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-nosave':
        saveplot = False
    elif arg == '-nocontour':
        plotcontours = False
    elif arg == '-noshow':
        showplot = False

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

# Get the time range in sec
t1 = translate_times(iter1, dirname, translate_from='iter')['val_sec']
t2 = translate_times(iter2, dirname, translate_from='iter')['val_sec']

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

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

# Compute the magnitude of B_m
bm = np.sqrt(br**2 + bt**2)

# Compute the streamfunction
psi = streamfunction(br, bt, rr, cost)

# Make CCW negative and CW positive
bm *= np.sign(psi)

if minmax is None:
    std = np.std(br)
    nstd = 3.
    minmax = -nstd*std, nstd*std
# Create plot
subplot_width_inches = 2.5
subplot_height_inches = 5.
margin_inches = 1/8
margin_top_inches = 1.25 # larger top margin to make room for titles

fig_width_inches = subplot_width_inches + 2*margin_inches
fig_height_inches = subplot_height_inches + margin_top_inches + margin_inches

fig_aspect = fig_height_inches/fig_width_inches
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
ax = fig.add_axes((margin_x, margin_y, subplot_width, subplot_height))

# Plot B_m magnitude
plot_azav (bm, rr, cost, fig=fig, ax=ax,\
    units = r'$\rm{G}$', plotcontours=False,\
    minmax=minmax)

# Plot streamfunction contours
#maxabs = np.max(np.abs(psi))
maxabs = 3*np.std(psi)
levels = np.linspace(-maxabs, maxabs, 20)
plot_azav (psi, rr, cost, fig=fig, ax=ax, plotfield=False, levels=levels)

# Label averaging interval
if rotation:
    time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
            + time_label + '\n' + (r'$\ (\Delta t = %.1f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'
else:
    time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.3f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'

# Make title
fsize = 12
line_height = 1./4./fig_height_inches
fig.text(margin_x, 1 - margin_y, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - margin_y - line_height,\
         r'$|\langle\mathbf{B}_m\rangle|$',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - margin_y - 2*line_height,\
        time_string, ha='left', va='top', fontsize=fsize, **csfont)

savefile = plotdir + dirname_stripped + '_Bm_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + '.png'
print ('Saving plot at %s ...' %savefile)
if saveplot:
    plt.savefig(savefile, dpi=300)
if showplot:
    plt.show()
