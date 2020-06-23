# Routine to compute signed and unsigned angular momentum in various
# radial zones
# By default the zones are whole CZ (for CZ-only runs)
# and CZ + RZ (for tachocline runs with domain_bounds specified)
# Created: 12/23/2018 (but really before)
#
# By default amom is calculated from the longest possible time
# interval in the AZ_Avgs data 
import numpy as np
import os, sys
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from get_parameter import get_parameter
from rayleigh_diagnostics import GridInfo
from get_eq import get_eq
from common import get_widest_range_file, strip_dirname, get_dict

# Get the name of the simulation directory
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)
datadir = dirname + '/data/'
print ("In directory ", dirname_stripped, "/")

# get density and grid info
eq = get_eq(dirname)
nr = eq.nr
rho = (eq.density).reshape((1, nr))

# get default AZ_Avgs file
the_file = get_widest_range_file(datadir, 'AZ_Avgs')

# By default integrate over the whole shell
rr = eq.radius

# By default, interate over different Chebyshev domains
try:
    domain_bounds = get_parameter(dirname, 'domain_bounds')
    ndom = len(domain_bounds) - 1
    r1 = domain_bounds[:ndom]
    r2 = domain_bounds[1:]
    ir1 = np.zeros(ndom, dtype='int')
    ir2 = np.zeros(ndom, dtype='int')
    for idom in range(ndom):
        ir1[idom] = np.argmin(np.abs(rr - r1[idom]))
        ir2[idom] = np.argmin(np.abs(rr - r2[idom]))
    # Remember the middle domain boundaries double up 
    # argmin gives first index of occurrence (higher radius)
    # This is OK for the ir1's, but the ir2's should be shifted to a
    # lower radius (up in index), except highest-radius (last) one
    ir2[:-1] += 1
    print ("domain_bounds in main_input detected, averaging over")
    print ("different Chebysehv domains by default")
except:
    r1 = np.array([get_parameter(dirname, 'rmin')])
    r2 = np.array([get_parameter(dirname, 'rmax')])
    ir1 = np.array([nr - 1])
    ir2 = np.array([0])
    ndom = 1
    print ("rmin/rmax in main_input detected, averaging over")
    print ("a single domain (probably the CZ)")

# Read in CLAs
args = sys.argv[2:]
nargs = len(args)
density = False

for i in range(nargs):
    arg = args[i]
    if arg == '-rmin':
        r1 = np.array(float(args[i+1]))
        ir1 = np.array([np.argmin(np.abs(rr - r1))])     
    elif arg == '-rmax':
        r2 = np.array(float(args[i+1]))
        ir2 = np.array([np.argmin(np.abs(rr - r2))])     
    elif arg == '-dens':
        density = True
    elif arg == '-usefile':
        the_file = args[i+1]
        the_file = the_file.split('/')[-1]

# Get the actual rmin/rmax values associated with the index
# For the default zones, this should have no effect
print(ir1)
r1 = rr[ir1]
r2 = rr[ir2]

# Get the AZ_Avgs data
print ('Getting data from ' + datadir + the_file)
di = get_dict(datadir + the_file)
vals = di['vals']
lut = di['lut']
xx = di['xx']
nt = di['nt']

# Get the azimuthal velocity
vp_av = vals[:, :, lut[3]]

# Get the radial/horizontal integration weights 
gi = GridInfo(dirname + '/grid_info', '')
rw = gi.rweights
tw = (gi.tweights).reshape((nt, 1))

# Get angular momentum (density) in the rotating frame. 
amom = rho*vp_av*xx
amom_u = np.abs(amom)

# Average over latitude
amom = np.sum(amom*tw, axis=0)
amom_u = np.sum(amom_u*tw, axis=0)

# Integrate the amom/unsigned amom over various domains
for idom in range(ndom):
    r1_loc = r1[idom]
    ir1_loc = ir1[idom]
    r2_loc = r2[idom]
    ir2_loc = ir2[idom]

    print("=====================================")
    print ("Zone #%i: %1.3e cm to %1.3e cm" %(idom + 1, r1_loc, r2_loc))
    #print("-------------------------------------")

    # Recall radius arrays are reversed, so ir2 < ir1
    amom_loc = np.sum((amom*rw)[ir2_loc:ir1_loc+1])/\
            np.sum(rw[ir2_loc:ir1_loc+1])
    amom_u_loc = np.sum((amom_u*rw)[ir2_loc:ir1_loc+1])/\
            np.sum(rw[ir2_loc:ir1_loc+1])
    if not density:
        shell_volume = 4./3.*np.pi*(r2_loc**3. - r1_loc**3.)
        amom_loc *= shell_volume
        amom_u_loc *= shell_volume

    if density:
        print ("%1.3e g/cm/s (total amom)" %amom_loc)
        print ("%1.3e g/cm/s (unsigned amom)" %amom_u_loc)
    else:
        print ("%1.3e g cm^2/s (total amom)" %amom_loc)
        print ("%1.3e g cm^2/s (unsigned amom)" %amom_u_loc)
print("=====================================")
