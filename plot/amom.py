# Author: Loren Matilsky
# Created: 06/11/2020
# This script generates angular momentum plotted in the meridional plane 
# for the Rayleigh run directory indicated by [dirname]. To use an AZ_Avgs file
# different than the one associated with the longest averaging range, use
# -usefile [complete name of desired AZ_Avgs file]
# Saves plot in
# [dirname]_amom_[first iter]_[last iter].png
# or 
# [dirname]_amom_spec_[first iter]_[last iter].png
# depending if -spec option was used (to plot specific amom, without density factor)

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
from get_eq import get_eq

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
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Set defaults
nlevs = 20
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
rbcz = None
minmax = None
rvals = None
spec = False

# Read in CLAs (if any) to change default variable ranges and other options
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-rbcz':
        rbcz = float(args[i+1])
    elif arg == '-nlevs':
        nlevs = int(args[i+1])
    elif arg == '-usefile':
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for rval_str in rvals_str:
            rvals.append(float(rval_str))
    elif arg == '-spec':
        spec = True
        
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
nr = di['nr']
nt = di['nt']

# Get the density
eq = get_eq(dirname)
rho = (eq.density).reshape((1, nr))

# Get angular momentum in the rotating frame. 
amom = vp_av*xx
if not spec:
    amom *= rho

# Create plot
subplot_width_inches = 2.5
subplot_height_inches = 5.
margin_inches = 1./8.
margin_top_inches = 1.25 # larger top margin to make room for titles
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
if spec:
    units = r'$\rm{cm^2\ s^{-1}}$'
else:
    units = r'$\rm{g\ cm^{-1}\ s^{-1}}$'
plot_azav (amom, rr, cost, fig=fig, ax=ax, units=units,\
        nlevs=nlevs, minmax=minmax, rvals=rvals)
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
if spec:
    amom_label = r'$\mathcal{L}\equiv r\sin\theta\langle v_\phi\rangle$'
else:
    amom_label = r'$\mathcal{L}\equiv \overline{\rho}r\sin\theta\langle v_\phi\rangle$'
fig.text(margin_x, 1 - margin_y - 2*line_height, amom_label,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - margin_y - 3*line_height, time_string,\
         ha='left', va='top', fontsize=fsize, **csfont)
if spec:
    basename = '_amom_spec_'
else:
    basename = '_amom_'
savefile = plotdir + dirname_stripped + basename + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + '.png'
print ('Saving plot at %s ...' %savefile)
plt.savefig(savefile, dpi=300)
plt.show()
