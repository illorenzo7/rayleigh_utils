import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append('/Users/loren/rayleigh/plot')
from azavg_util import plot_azav, streamfunction
from binormalized_cbar import MidpointNormalize
from diagnostic_reading import ReferenceState
from get_parameter import get_parameter
from subprocess import call
from common import get_widest_range_file, get_iters_from_file

# Read in the name of the run directory
dirname = sys.argv[1]
dirname_stripped = dirname.split('/')[-1]
#dirname = '/castor/loma3853/rayleigh/n3_3/'

# Get data and plot directories
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Read in eflux_radial data file
eflux_file = get_widest_range_file(datadir, 'eflux_radial_merplane')
vflux, kflux, cflux, hflux, eflux, tflux = np.load(datadir + eflux_file)

# Get grid info
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr, nt = len(rr), len(tt)

# Read command-line arguments (CLAs)
my_boundstype = 'minmax'
user_specified_minmax = False
showplot = False
my_nlevs = 20
my_min, my_max = -10, 10 # placeholders

var_to_plot = cflux
varname = 'cflux'

titles = {'cflux': 'Conductive Flux', 'vflux': 'Viscous Flux',\
          'kflux': 'Kinetic Energy Flux', 'eflux': 'Enthalpy Flux'}

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-minmax'):
        my_boundstype = 'manual'
        my_min, my_max = float(args[i+1]), float(args[i+2])
        user_specified_minmax = True
    elif (arg == '-show'):
        showplot = True
    elif (arg == '-nlevs'):
        my_nlevs = int(args[i+1])
    elif (arg == '-visc'):
        var_to_plot = vflux
        varname = 'vflux'
    elif (arg == '-ke'):
        var_to_plot = kflux
        varname = 'kflux'
    elif (arg == '-enth'):
        var_to_plot = eflux
        varname = 'eflux'
        
fig, ax = plt.subplots()
plot_azav (fig, ax, var_to_plot, rr, cost, sint,\
           units = r'$\rm{ergs}/\rm{s}/\rm{cm}^2$', plotcontours=False,\
           boundstype = my_boundstype, caller_minmax = (my_min, my_max))

plt.title (titles[varname], fontsize=16)
plt.tight_layout()

iter1, iter2 = get_iters_from_file(eflux_file)

savename = dirname_stripped + '_' + varname + '_merplane_' + str(iter1).zfill(8) + '_' +\
    str(iter2).zfill(8) + '.png'
savefile = plotdir + savename

print ('Saving file at ' + savefile + ' ...')

plt.savefig(savefile, dpi=300)
#plt.close()

if (showplot):
#    call(['open', savefile])   
    plt.show()