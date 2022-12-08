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
from common import *

# Get the name of the run directory
dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist
datadir = dirname + '/data/'
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

radatadir = dirname + '/Meridional_Slices/'

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Read in CLAs
args = sys.argv[2:]
nargs = len(args)

if (nargs == 0):
    index_first, index_last = nfiles - 101, nfiles - 1  
    # By default average over the last 100 files
else:
    index_first, index_last = get_desired_range(int_file_list, args)

# Set the timeavg savename by the directory, what we are saving, and first and last
# iteration files for the average
savename = dirname_stripped + '_enstrophy_from_mer_' +\
        file_list[index_first] + '_' + file_list[index_last] + '.pkl'
savefile = datadir + savename    

# Get grid info from first mer slice file
mer0 = Meridional_Slices(radatadir + file_list[index_first], '')
rr = mer0.radius
ri, ro = np.min(rr), np.max(rr)
d = ro - ri
rr_depth = (ro - rr)/d
rr_height = (rr - ri)/d
sint = mer0.sintheta
cost = mer0.costheta
tt = np.arccos(cost)
tt_lat = (np.pi/2 - tt)*180/np.pi
nr = mer0.nr
nt = mer0.ntheta
phivals = mer0.phi
# phiinds = mer0.phi_indices 
# This attribute of Meridional_Slices does not actually seem to be there!
nphi = mer0.nphi
phivals_lon = 180*phivals/np.pi
# compute some derivative quantities for the grid
tt_2d, rr_2d = np.meshgrid(tt, rr, indexing='ij')
sint_2d = np.sin(tt_2d); cost_2d = np.cos(tt_2d)
xx = rr_2d*sint_2d
zz = rr_2d*cost_2d

# Initialize the count (will rise as the averaging commences)
count = 0
iter1, iter2 = int_file_list[index_first], int_file_list[index_last]

# Initialize empty "vals" array for the time average
vals = np.zeros((nphi, nt, nr, 3))

# Average over the relevant data range, summing everything and then dividing
#   by the number of "slices" added at the end
print ('Considering Meridional_Slices files %s through %s for the average ...'\
       %(file_list[index_first], file_list[index_last]))
for i in range(index_first, index_last + 1):
    print ('Adding Meridional_Slices/%s to the average ...' %file_list[i])
    if i == index_first:
        mer = mer0
    else:   
        mer = Meridional_Slices(radatadir + file_list[i], '')

    local_ntimes = mer.niter
    for j in range(local_ntimes):
        vals += (mer.vals[:, :, :,\
                [mer.lut[301], mer.lut[302], mer.lut[303]], j])**2
        count += 1

vals /= count
print ('Averaged over %i Meridional_Slices ...' %count)

# Save the avarage
print ('Saving file at ' + savefile + ' ...')
f = open(savefile, 'wb')
pickle.dump({'vals': vals, 'lut': mer0.lut, 'count': count, 'iter1': iter1, 'iter2': iter2, 'phivals': phivals, 'phivals_lon': phivals_lon, 'nphi': nphi, 'rr': rr, 'rr_depth': rr_depth, 'rr_height': rr_height, 'nr': nr, 'ri': ri, 'ro': ro, 'd': d, 'tt': tt, 'tt_lat': tt_lat, 'sint': sint, 'cost': cost,'nt': nt, 'rr_2d': rr_2d, 'tt_2d': tt_2d, 'sint_2d': sint_2d, 'cost_2d': cost_2d, 'xx': xx, 'zz': zz}, f, protocol=4)
f.close()
