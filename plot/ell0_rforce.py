###############################################
# Author: Loren Matilsky
# Date created: 02/19/2018
#
# This script plots the radial forces (ell=0 components) 
# as functions of radius using from the Shell_Avgs data
# of custom outputs 2216--2220 (2223 if magnetism)

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
from common import get_widest_range_file, strip_dirname, get_iters_from_file
from get_parameter import get_parameter

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Find the Shell_Avgs file(s) in the data directory. If there are multiple, by
# default choose the one with widest range in the average
Shell_Avgs_file = get_widest_range_file(datadir, 'Shell_Avgs')

# Get command-line arguments to adjust the interval of averaging files
user_specified_minmax = False
user_specified_rnorm = False
try:
    magnetism = get_parameter(dirname, 'magnetism')
except:
    magnetism = False
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-usefile'):
        Shell_Avgs_file = args[i+1]
        Shell_Avgs_file = Shell_Avgs_file.split('/')[-1]
    elif (arg == '-minmax'):
        user_specified_minmax = True
        my_min, my_max = float(args[i+1]), float(args[i+2])
    elif (arg == '-rnorm'):
        user_specified_rnorm = True
        user_supplied_rnorm = float(args[i+1])

#Create the plot
lw = 1. # regular-width lines

# Read in the flux data
print ('Getting ell=0 radial forces from ' + datadir + Shell_Avgs_file + ' ...')
try:
    di = np.load(datadir + Shell_Avgs_file, encoding='latin1').item()
except:
    f = open(Shell_Avgs_file, 'rb')
    di = pickle.load(f)
    f.close()

vals = di['vals']
lut = di['lut']
iter1, iter2 = di['iter1'], di['iter2']
rr = di['rr']

# Make the plot name, labelling the first/last iterations we average over
savename = dirname_stripped + '_ell0_rforce_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

qindex_buoy = lut[1246]
qindex_prs = lut[1247]
qindex_visc = lut[2216]
qindex_cor = lut[2217]
qindex_adv_tot = lut[2218]
qindex_adv_fluc = lut[2219]
qindex_adv_mean = lut[2220]

if magnetism:
    qindex_lor_tot = lut[2221]
    qindex_lor_fluc = lut[2222]
    qindex_lor_mean = lut[2223]

buoy = vals[:, qindex_buoy]
prs = vals[:, qindex_prs]
visc = vals[:, qindex_visc]
cor = vals[:, qindex_cor]
adv_tot = vals[:, qindex_adv_tot]
adv_fluc = vals[:, qindex_adv_fluc]
adv_mean = vals[:, qindex_adv_mean]
tforce = buoy + prs + visc + cor + adv_tot

if magnetism:
    lor_tot = vals[:, qindex_lor_tot]
    lor_fluc = vals[:, qindex_lor_fluc]
    lor_mean = vals[:, qindex_lor_mean]
    tforce += lor_tot
    
# Create the plot
Rsun = 6.955e10 # Normalize the radius by the solar radius

# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if not user_specified_rnorm:
    rr_n = rr/Rsun
else:
    rr_n = rr/user_supplied_rnorm                                           

plt.plot(rr_n, buoy, label=r'$(\mathbf{f}_{\rm{buoy}})_{l=0}$', linewidth=lw)
plt.plot(rr_n, prs, label=r'$(\mathbf{f}_{\rm{prs}})_{l=0}$', linewidth=lw)
plt.plot(rr_n, visc, label=r'$(\mathbf{f}_{\rm{visc}})_{l=0}$', linewidth=lw)
plt.plot(rr_n, cor, label=r'$(\mathbf{f}_{\rm{cor}})_{l=0}$', linewidth=lw)
plt.plot(rr_n, adv_tot, 'k', label=r'$(\mathbf{f}_{\rm{adv}})_{l=0}$', linewidth=lw)
plt.plot(rr_n, adv_fluc, label=r'$(\mathbf{f}^\prime_{\rm{adv}})_{l=0}$', linewidth=lw)
plt.plot(rr_n, adv_mean, label=r'$(\overline{\mathbf{f}}_{\rm{adv}})_{l=0}$', linewidth=lw)
plt.plot(rr_n, adv_mean + adv_fluc, 'r--', \
label=r'$(\overline{\mathbf{f}}_{\rm{adv}})_{l=0} + (\mathbf{f}^\prime_{\rm{adv}})_{l=0}$',\
plt.plot(rr_n, tforce, label=r'$(\mathbf{f}_{\rm{tot}})_{l=0}$', linewidth=lw)
linewidth=lw)
if magnetism:
    plt.plot(rr_n, lor_tot, label=r'$(\mathbf{f}_{\rm{lor}})_{l=0}$', linewidth=lw)
    plt.plot(rr_n, lor_fluc, label=r'$(\mathbf{f}^\prime_{\rm{lor}})_{l=0}$', linewidth=lw)
    plt.plot(rr_n, lor_mean, label=r'$(\overline{\mathbf{f}}_{\rm{lor}})_{l=0}$', linewidth=lw)
    plt.plot(rr_n, lor_mean + lor_fluc, 'b--', \
    label=r'$(\overline{\mathbf{f}}_{\rm{lor}})_{l=0} + (\mathbf{f}^\prime_{\rm{lor}})_{l=0}$',\
    linewidth=lw)

# Get the y-axis in scientific notation
plt.ticklabel_format(useMathText=True, axis='y', scilimits=(0,0))

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Set the x limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
delta_x = xmax - xmin
plt.xlim(xmin, xmax)

# Set the y-limits (the following values seem to "work well" for my models
# so far...perhaps adjust this in the future. 
#ymin, ymax = -0.7, 1.3
#if user_specified_minmax:
#    ymin, ymax = my_min, my_max
#delta_y = ymax - ymin
#plt.ylim(ymin, ymax)

# Label the axes
if not user_specified_rnorm:
    plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
else:
    plt.xlabel(r'r/(%.1e cm)' %user_supplied_rnorm, fontsize=12, **csfont)
plt.ylabel('force density',\
        fontsize=12, **csfont)

# Make title
plt.title(dirname_stripped + '\n' + 'l=0 radial force, ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8), **csfont)

# Create a see-through legend
plt.legend(loc='lower left', shadow=True, ncol=3, fontsize=10)

# Last command
plt.tight_layout()

# Save the plot
print ('Saving the ell=0 rforce plot at ' + plotdir + savename + ' ...')
plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()
