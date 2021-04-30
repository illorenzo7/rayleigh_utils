# Author: Loren Matilsky
# Created: 10/31/2019
# This script generates the spherically averaged convective Reynolds 
# number (Re), plotted along radial lines for
# the Rayleigh run directory indicated by [dirname]. 
# Computes contributions from v_r, v_theta, and v_phi separately 
# To use  time-averaged 
# AZ_Avgs file different than the one associated with the longest averaging 
# range, use -usefile [complete name of desired vavg file]
# Saves plot in
# [dirname]_Re_rslice_[first iter]_[last iter].png

# Import relevant modules
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import *

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# get rho
eq = get_eq(dirname)
rho = eq.density

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'

# Set defaults
rnorm = None
minmax = None
logscale = False
rvals = [] # user can specify radii to mark by vertical lines
tag = ''
use_hrho = False
the_file = get_widest_range_file(datadir, 'Shell_Avgs')

# Read command-line arguments (CLAs)
plotdir = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-usefile':
        the_file = args[i+1]
        the_file = the_file.split('/')[-1]
    elif arg == '-rnorm':
        rnorm = float(args[i+1])
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-log':
        logscale = True
    elif arg == '-hrho':
        use_hrho = True
    elif arg == '-tag':
        tag = '_' + args[i+1]
    elif arg == '-depths':
        strings = args[i+1].split()
        for st in strings:
            rval = ro - float(st)*d
            rvals.append(rval)
    elif arg == '-depthscz':
        rm = domain_bounds[1]
        dcz = ro - rm
        strings = args[i+1].split()
        for st in strings:
            rval = ro - float(st)*dcz
            rvals.append(rval)
    elif arg == '-depthsrz':
        rm = domain_bounds[1]
        drz = rm - ri
        strings = args[i+1].split()
        for st in strings:
            rval = rm - float(st)*drz
            rvals.append(rval)
    elif arg == '-rvals':
        rvals = []
        strings = args[i+1].split()
        for st in strings:
            rval = float(st)*rsun
            rvals.append(rval)
    elif arg == '-rvalscm':
        rvals = []
        strings = args[i+1].split()
        for st in strings:
            rval = float(st)
            rvals.append(rval)

# Read in vavg data
print ('Reading Shell_Avgs data from ' + the_file + ' ...')
di = get_dict(the_file)
vals = di['vals']
lut = di['lut']
iter1, iter2 = get_iters_from_file(the_file)

di_grid = get_grid_info(dirname)
rr = di_grid['rr']

# Derivative grid info
nr = len(rr)
ri, ro = np.min(rr), np.max(rr)
shell_depth = ro - ri

# Convective B amplitudes, get these from ME
rme = vals[:, 0, lut[1102]]
tme = vals[:, 0, lut[1103]]
pme = vals[:, 0, lut[1104]]

frme = vals[:, 0, lut[1110]]
ftme = vals[:, 0, lut[1111]]
fpme = vals[:, 0, lut[1112]]

mrme = rme - frme
mtme = tme - ftme
mpme = pme - fpme

eightpi = 8.*np.pi

bsq_r = rme*eightpi
bsq_t = tme*eightpi
bsq_p = pme*eightpi
bsq = bsq_r + bsq_t + bsq_p

fbsq_r = frme*eightpi
fbsq_t = ftme*eightpi
fbsq_p = fpme*eightpi
fbsq = fbsq_r + fbsq_t + fbsq_p

mbsq_r = mrme*eightpi
mbsq_t = mtme*eightpi
mbsq_p = mpme*eightpi
mbsq = mbsq_r + mbsq_t + mbsq_p

amp_b = np.sqrt(bsq)
amp_br = np.sqrt(bsq_r)
amp_bt = np.sqrt(bsq_t)
amp_bp = np.sqrt(bsq_p)

famp_b = np.sqrt(fbsq)
famp_br = np.sqrt(fbsq_r)
famp_bt = np.sqrt(fbsq_t)
famp_bp = np.sqrt(fbsq_p)

mamp_b = np.sqrt(mbsq)
mamp_br = np.sqrt(mbsq_r)
mamp_bt = np.sqrt(mbsq_t)
mamp_bp = np.sqrt(mbsq_p)

# Create the plot
fig = plt.figure()
ax = fig.add_subplot(111)

# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if rnorm is None:
    rr_n = rr/rsun
else:
    rr_n = rr/rnorm                                           

# Plot B amps vs radius
colors = ['k', 'r', 'g', 'b', 'm', 'c']
ax.plot(rr_n, amp_br/amp_bt, colors[0] + '-', label='tot')

ax.plot(rr_n, famp_br/famp_bt, colors[1] + '--', label='m != 0')

ax.plot(rr_n, mamp_br/mamp_bt, colors[2] + ':', label='m = 0')

# Label the axes
if rnorm is None:
    plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
else:
    plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12, **csfont)

plt.ylabel('rms b_r / rms b_t',fontsize=12,\
        **csfont)

# Set the axis limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
plt.xlim((xmin, xmax))
xvals = np.linspace(xmin, xmax, 100)

if not minmax is None:
    ymin, ymax = minmax
    plt.ylim((ymin, ymax))

if plotdir is None:
    plotdir = dirname + '/plots/'
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)

# Create a title    
plt.title(dirname_stripped + '\n' +'B Field Amplitudes, ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8), **csfont)
plt.legend()

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')
plt.tight_layout()

savefile = plotdir + dirname_stripped + '_br_to_bt_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'
print('Saving plot at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.show()
