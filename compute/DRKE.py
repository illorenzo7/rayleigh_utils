# Routine to figure out how much energy is in the equilibrated 
# differential rotation, normalized or not by the energy in the frame
# rotation
# Created: 10/20/2020

import numpy as np
import os, sys
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from get_parameter import get_parameter
from get_eq import get_eq
from time_scales import compute_Prot, compute_tdt
from rayleigh_diagnostics import GridInfo
from common import get_file_lists, get_widest_range_file, strip_dirname,\
        get_dict

# Main directory and data directory
dirname = sys.argv[1]
datadir = dirname + '/data/'

# Read in CLAs
args = sys.argv[2:]
nargs = len(args)

# Find the etrace file(s) in the data directory. If there are multiple, by
# default choose the one with widest range in the trace.
the_file = get_widest_range_file(datadir, 'G_Avgs')

# Set defaults
xiter = False
from0 = False
magnetism = False
ylog = False
minmax = None
xminmax = None
xmin = None
xmax = None
plot_inte = False
plot_tote = False
savename = None
tol = 0.05
plot_equil_time = False
chosen_eqtime = None # choose eq time a priori

# Get command-line arguments
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-xiter': # plot w.r.t. iterations
        xiter = True
    elif arg == '-gtr':
        the_file = get_widest_range_file(datadir, 'trace_G_Avgs')
        iterstart_desired = float(args[i+1])
    elif arg == '-usefile':
        the_file = args[i+1]
        the_file = the_file.split('/')[-1]
    elif arg == '-from0':
        from0 = True
    elif arg == '-mag':
        magnetism = True
    elif arg == '-log':
        ylog = True
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-xmin':
        xmin = float(args[i+1])
    elif arg == '-xmax':
        xmax = float(args[i+1])
    elif arg == '-inte':
        plot_inte = True
    elif arg == '-tote':
        plot_tote = True
    elif arg == '-name':
        savename = args[i+1] + '.png'
    elif arg == '-tol':
        tol = float(args[i+1])
    elif arg == '-equil':
        plot_equil_time = True
    elif arg == '-eqtime':
        chosen_eqtime = float(args[i+1])


# get density and grid info
eq = get_eq(dirname)
rho = eq.rho
gi = GridInfo(dirname + '/grid_info', '')
rr = gi.radius
rw = gi.rweights
nr = len(rr)

# Get the frame rate energy
Om0 = 2*np.pi/compute_Prot(dirname)

# Get the radial integration weights and density
rw = gi.rweights

w0 = Om0**2.0*np.sum(rho*rr**2.0*rw)/3.0

# Compute the equilibrated DRKE
# Read in the KE data (dictionary form)
if 'trace' in the_file:
    print ('Getting trace of KEs from ' + datadir + the_file)
else:
    print ('Getting average KEs from ' + datadir + the_file)
di = get_dict(datadir + the_file)
vals = di['vals']
lut = di['lut']
if 'trace' in the_file:
    iters = di['iters']
    istart = np.argmin(np.abs(iters - iterstart_desired))
    iterstart = iters[istart]
    pke = np.mean(vals[lut[404]][istart:])
    fpke = np.mean(vals[lut[412]][istart:])
    print ("Averaging trace from iters")
    print("%08i to %08i" %(iterstart, iters[-1]))
else:
    pke = vals[lut[404]]
    fpke = vals[lut[412]]
mpke = pke - fpke
DRKE = mpke

print ("frame KE = %1.2e" %w0)
print ("DRKE = %1.2e" %DRKE)
print ("ratio = %1.2e" %(DRKE/w0))
