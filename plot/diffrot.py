# Author: Loren Matilsky
# Created: 05/14/2018
# This script generates differential rotation plotted in the meridional plane 
# for the Rayleigh run directory indicated by [dirname]. To use an AZ_Avgs file
# different than the one associated with the longest averaging range, use
# -usefile [complete name of desired AZ_Avgs file]
# Saves plot in
# [dirname]_diffrot_[first iter]_[last iter].npy

import numpy as np
import pickle
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from azav_util import plot_azav
from common import get_widest_range_file, strip_dirname, get_dict
from get_parameter import get_parameter
from time_scales import compute_Prot, compute_tdt
from translate_times import translate_times

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Split dirname_stripped into two lines if it is very long
if len(dirname_stripped) > 25:
    dirname_stripped_title = dirname_stripped[:25] + '\n' +\
            dirname_stripped[25:]
else:
    dirname_stripped_title = dirname_stripped

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Set defaults
my_boundstype = 'manual'
user_specified_minmax = False 
my_nlevs = 20
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
rbcz = None

# Read in CLAs (if any) to change default variable ranges and other options
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-minmax'):
        my_min, my_max = float(args[i+1]), float(args[i+2])
        user_specified_minmax = True
    elif arg == '-rbcz':
        rbcz = float(args[i+1])
    elif (arg == '-nlevs'):
        my_nlevs = int(args[i+1])
    elif (arg == '-usefile'):
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
        
# Read in AZ_Avgs data
print ('Getting data from ' + datadir + AZ_Avgs_file + ' ...')
di = get_dict(datadir + AZ_Avgs_file)

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']

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

vr_av, vt_av, vp_av = vals[:, :, lut[1]], vals[:, :, lut[2]],\
        vals[:, :, lut[3]]

# Get necessary grid info
rr = di['rr']
cost = di['cost']
sint = di['sint']
tt_lat = di['tt_lat']
xx = di['xx']

# Get differential rotation in the rotating frame. 
Om = vp_av/xx
diffrot = Om*1.0e9/2/np.pi # rad/s --> nHz

# Maximum differential rotation over whole meridional plane
it15, it75 = np.argmin(np.abs(tt_lat - 15)), np.argmin(np.abs(tt_lat - 75))
     # ignore problematic poles 
global_min, global_max = np.min(diffrot[it15:it75, :]),\
     np.max(diffrot[it15:it75, :])
Delta_Om = global_max - global_min
maxabs = np.max((np.abs(global_min), np.abs(global_max)))

if (not user_specified_minmax):
    my_min, my_max = -maxabs, maxabs

# Create plot
subplot_width_inches = 2.5
subplot_height_inches = 5.
margin_inches = 1./8.
margin_top_inches = 2. # larger top margin to make room for titles
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)

fig_width_inches = subplot_width_inches + 2*margin_inches
fig_height_inches = subplot_height_inches + margin_top_inches +\
        margin_bottom_inches

fig_aspect = fig_height_inches/fig_width_inches
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
ax = fig.add_axes((margin_x, margin_bottom, subplot_width, subplot_height))

plot_azav (diffrot, rr, cost, fig=fig, ax=ax, units='nHz',\
        nlevs=my_nlevs, minmax = (my_min, my_max))

# Make title + label diff. rot. contrast and no. contours
# Label averaging interval
if rotation:
    time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
            + time_label + '\n' + (r'$\ (\Delta t = %.1f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'
else:
    time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.3f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'
fsize = 12.
line_height = 1./4./fig_height_inches
fig.text(margin_x, 1 - margin_y, dirname_stripped_title,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - margin_y - 2*line_height, r'$\Omega - \Omega_0$',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - margin_y - 3*line_height, time_string,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - margin_y - 5*line_height,\
         r'$\Delta\Omega_{\rm{tot}} = %.1f\ nHz$' %Delta_Om,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - margin_y - 6*line_height,\
         'nlevs = %i' %my_nlevs,
         ha='left', va='top', fontsize=fsize, **csfont)
savefile = plotdir + dirname_stripped + '_diffrot_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + '.png'
print ('Saving plot at %s ...' %savefile)
plt.savefig(savefile, dpi=300)
plt.show()
