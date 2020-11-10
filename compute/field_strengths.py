# Routine to average Rayleigh Meridional_Slices ENSTROPHY data in time
# Created by: Loren Matilsky
# On: 04/10/2019
##################################################################
# This routine computes the average in time of the ENSTROPHY values in the
# Meridional_Slices data for a particular simulation. 

# By default, the routine averages over the last 100 files of datadir, though
# the user can specify a different range in sevaral ways:
# -n 10 (last 10 files)
# -range iter1 iter2 (names of start and stop data files; 
# names can also be "first" or "last")
# -centerrange iter0 nfiles (average about central file iter0 over nfiles)

# Import relevant modules
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import Meridional_Slices
from common import get_file_lists, get_desired_range

# Get the name of the run directory
dirname = sys.argv[1]

radatadir = dirname + '/Meridional_Slices/'

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Read in CLAs
args = sys.argv[2:]
nargs = len(args)

the_tuple = get_desired_range(int_file_list, args)
if the_tuple is None: # By default plot the last 10 Shell_Slices
    index_first, index_last = 0, nfiles - 1  
    # By default search strenghts over over the last 100 files
else:
    index_first, index_last = the_tuple
nused = index_last - index_first + 1

# Defaults
rbcz = None # by default treat strengths as a whole CZ
        # (if repeated radii are detected, will change rbcz
        # to the repeated radius (if not already specified))
nskip = nused//11
# by default look at 11 meridional slice

for i in range(nargs):
    arg = args[i]
    if arg == '-ntot': 
        ntot = int(args[i+1])
        nskip = nused//ntot
    elif arg == '-rbcz':
        rbcz = float(args[i+1])

# Get grid info from first mer slice file
mer0 = Meridional_Slices(radatadir + file_list[index_last], '')
rr = mer0.radius
nr = mer0.nr
# Search for repeated radii
tol = 1.0e-12
for ir in range(1, nr):
    if np.abs(rr[ir] - rr[ir-1]) < tol:
        if rbcz is None:
            rbcz = rr[ir - 1]
        print ("Repeated radius detected: r = %8.2e" %rr[ir-1])
if not rbcz is None:
    irbcz = np.argmin(np.abs(rr - rbcz))
    nrcz = irbcz + 1
    nrrz = nr - nrcz

sint = mer0.sintheta
cost = mer0.costheta
tt = np.arccos(cost)
tt_lat = (np.pi/2 - tt)*180/np.pi
nt = mer0.ntheta
nphi = mer0.nphi
#phivals = mer0.phi
# compute some derivative quantities for the grid
#tt_2d, rr_2d = np.meshgrid(tt, rr, indexing='ij')
#sint_2d = np.sin(tt_2d); cost_2d = np.cos(tt_2d)
#xx = rr_2d*sint_2d
#zz = rr_2d*cost_2d

# Search for B-field extrema over the relevant data range,
print ("Considering Meridional_Slices files %s through %s for max/min B"\
       %(file_list[index_first], file_list[index_last]))
if not rbcz is None:
    print ("nr = ", nr)
    print ("using rbcz = %7.2e, ir = %i" %(rbcz, irbcz))
    print ("nrcz = ", nrcz)
    print ("nrrz = ", nrrz)

# Initialize values to store maxima
if rbcz is None:
    bp = -np.inf
    bpm = -np.inf
    bpp = -np.inf
else:
    bprz = -np.inf
    bpmrz = -np.inf
    bpprz = -np.inf
    bpcz = -np.inf
    bpmcz = -np.inf
    bppcz = -np.inf

for i in range(index_last, index_first - 1, -nskip):
    print ('Scanning Meridional_Slices/%s for maxes/mins' %file_list[i])
    if i == index_last:
        mer = mer0
    else:   
        mer = Meridional_Slices(radatadir + file_list[i], '')

    local_ntimes = mer.niter
    for j in range(local_ntimes):
        if rbcz is None:
            bp_loc = mer.vals[:, :, :, mer.lut[803], j]
            bpm_loc = np.mean(bp_loc, axis=0).reshape((1,nt,nr))
            bpp_loc = bp_loc - bp_loc
            bp = max(bp, np.max(np.abs(bp_loc)))
            bpm = max(bpm, np.max(np.abs(bpm_loc)))
            bpp = max(bpp, np.max(np.abs(bpp_loc)))
        else:
            bprz_loc = mer.vals[:, :, irbcz + 1:, mer.lut[803], j]
            bpmrz_loc = np.mean(bprz_loc, axis=0).reshape((1, nt, nrrz))
            bpprz_loc = bprz_loc - bpmrz_loc

            bpcz_loc = mer.vals[:, :, :irbcz, mer.lut[803], j]
            bpmcz_loc = np.mean(bpcz_loc, axis=0).reshape((1, nt, nrcz))
            bppcz_loc = bpcz_loc - bpmcz_loc

            bprz = max(bprz, np.max(np.abs(bprz_loc)))
            bpcz = max(bpcz, np.max(np.abs(bpcz_loc)))

            bpmrz = max(bpmrz, np.max(np.abs(bpmrz_loc)))
            bpmcz = max(bpmcz, np.max(np.abs(bpmcz_loc)))

            bpprz = max(bpprz, np.max(np.abs(bpprz_loc)))
            bppcz = max(bppcz, np.max(np.abs(bppcz_loc)))

fmt = "%7.2e"
print ("-------------------------------------------------")
if rbcz is None:
    print (("max |B_phi|        = " + fmt + " G") %bp)
    print (("max |<B_phi>|      = " + fmt + " G") %bpm)
    print (("max |B_phi'|       = " + fmt + " G") %bpp)
else:
    print (("max |B_phi| (RZ)   = " + fmt + " G") %bprz)
    print (("max |<B_phi>| (RZ) = " + fmt + " G") %bpmrz)
    print (("max |B_phi'| (RZ)  = " + fmt + " G") %bpmrz)
    print ("-------------------------------------------------")
    print (("max |B_phi|   (CZ) = " + fmt + " G") %bpcz)
    print (("max |<B_phi>| (RZ) = " + fmt + " G") %bpmrz)
    print (("max |B_phi'|  (RZ)  = " + fmt + " G") %bpmrz)
print ("-------------------------------------------------")
