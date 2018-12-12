import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
from rayleigh_diagnostics import AZ_Avgs
from common import get_file_lists
dirname = sys.argv[1]
datadir = dirname + '/data/'

# Create 'datadir' if it doesn't exist already
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

# Get all the file names in datadir and their integer counterparts
azdir = dirname + '/AZ_Avgs/'
radatadir = dirname + '/AZ_Avgs/'
file_list, int_file_list, nfiles = get_file_lists(radatadir)

file_to_use = file_list[-1]
nargs = len(sys.argv[2:])
args = sys.argv[2:]
for i in range(nargs):
    arg = args[i]
    if arg == '-use': # use a different iteration number
        iter_to_use = int(args[i+1])
        file_to_use = file_list[np.argmin(np.abs(int_file_list - iter_to_use))]

print ('Getting grid_info from AZ_Avgs/%s ...' %file_to_use)
az0 = AZ_Avgs(azdir + file_to_use, '')

rr = az0.radius
cost = az0.costheta
sint = az0.sintheta
tt = np.arccos(cost)

# Shell dimensions
ro = rr[0]
ri = rr[len(rr) - 1]
d = ro - ri
# radial grid in "percent depth"
rr_depth = (ro - rr)/d

np.save(datadir + 'grid_info.npy',(rr,tt,cost,sint,rr_depth,ri,ro,d))

nr, nt = len(rr), len(tt)
print ('nr = %i, nth = %i' %(nr, nt))
