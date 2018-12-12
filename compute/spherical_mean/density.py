import numpy as np
import sys
import os
from diagnostic_reading import ShellAverage, ReferenceState
from get_parameter import get_parameter

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

pressure = np.zeros(nr)
entropy = np.zeros(nr)

count = 0
navg = 100
for ii in range(nfiles-navg, nfiles):
    print ('Adding G_Avgs/%s to the average...' %files[ii])
    shellavg = ShellAverage(radatadir + files[ii], '')
    local_ntimes = shellavg.niter
    for jj in range(local_ntimes):
        pressure_loc = shellavg.vals[:,0,shellavg.lut[65],jj]
        entropy_loc = shellavg.vals[:,0,shellavg.lut[64], jj]
        pressure += pressure_loc
        entropy += entropy_loc
        count += 1

pressure /= count
entropy /= count

# Get reference values
cp = get_parameter(dirname, 'pressure_specific_heat')
poly_n = get_parameter(dirname, 'poly_n')
gamma = 1. + 1./poly_n
ref = ReferenceState(dirname + '/reference', '')
p_ref = ref.pressure
rho_ref = ref.density

density_ratio = pressure/p_ref/gamma - entropy/cp
density = density_ratio*rho_ref

savefile = datadir + 'rho_spherical_mean.npy'
np.save(savefile, density)
