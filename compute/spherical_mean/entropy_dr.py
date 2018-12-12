import numpy as np
import sys
import os
from diagnostic_reading import ShellAverage

dirname = sys.argv[1]
radatadir = dirname + '/Shell_Avgs/'

datadir = dirname + '/data/'
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)


files = os.listdir(radatadir)
nfiles = len(files)
files.sort()

shellavg0 = ShellAverage(radatadir + files[0], '')
nr = shellavg0.nr

entropy_dr = np.zeros(nr)

count = 0
navg = 100
for ii in range(nfiles-navg, nfiles):
    print ('Adding Shell_Avgs/%s to the average...' %files[ii])
    shellavg = ShellAverage(radatadir + files[ii], '')
    local_ntimes = shellavg.niter
    for jj in range(local_ntimes):
        entropy_dr_loc = shellavg.vals[:,0,shellavg.lut[70],jj]
        entropy_dr += entropy_dr_loc
        count += 1

entropy_dr /= count

savefile = datadir + 's_dr_spherical_mean.npy'
np.save(savefile, entropy_dr)
