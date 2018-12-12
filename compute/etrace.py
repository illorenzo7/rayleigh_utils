# Author: Loren Matilsky
# Date created: 02/09/2017
#
#  Plots different KE averages for each recorded timestep
#
#  This example routine makes use of the GlobalAverage
#  data structure associated with the G_Avg output.
#  Upon initializing a GlobalAverage object, the 
#  object will contain the following attributes:
#
#    ----------------------------------
#    self.nrec                  : number of time steps
#    self.nq                     : number of diagnostic quantities output
#    self.qv[0:nq-1]             : quantity codes for the diagnostics output
#    self.vals[0:nrec-1,0:nq-1] : Globally averaged diagnostics as function of time and quantity index
#    self.iters[0:nrec-1]       : The time step numbers stored in this output file
#    self.time[0:nrec-1]        : The simulation time corresponding to each time step
#    self.lut                    : Lookup table for the different diagnostics output

# Import relevant modules
import numpy as np
import sys, os
from diagnostic_reading import GlobalAverage
from common import get_file_lists, get_desired_range, strip_dirname

# Get the name of the run directory
dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist
datadir = dirname + '/data/'
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

# Directory where the Rayleigh data is kept
data_type = 'G_Avgs'
radatadir = dirname + '/' + data_type + '/'

file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Get grid info (if it's not computed already using grid_info.py, this will fail)
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr = len(rr)
nt = len(tt)

# Read in CLAs
args = sys.argv[2:]
nargs = len(args)

if (nargs == 0):
    index_first, index_last = 0, nfiles   # By default do an e-trace over the 
                                          # whole run
else:
    index_first, index_last = get_desired_range(int_file_list, args)

# Set the savename by the directory, what we are saving, and first and last
# iteration files for the average
savename = dirname_stripped + '_etrace_' + file_list[index_first] + '_' +\
    file_list[index_last - 1] + '.npy'
savefile = datadir + savename    

#Initialize some empy lists to hold various KE variables

# KE, broken up into total, radial, latitudinal, azimuthal
ke = []
rke = []
tke = []
pke = []

# KE due to mean motion (DR + MC)
mke = []
mrke = []
mtke = []
mpke = []

# KE due to fluctuations (convection)
fke = []
frke = []
ftke = []
fpke = []

# We also define an empty list to hold the time. This will be the time
# since the beginning of the simulation (first output file).

times = [] 
iters = []

# Next, we loop over all files and grab the desired data
# from the file, storing it in appropriate lists as we go.

#for i in range(10): # for debugging purposes
for i in range(index_first, index_last):
    # read in the files needed one by one
    a = GlobalAverage(file_list[i],path=radatadir)
    print ('Adding %s to the time series ...' %(data_type + '/' + file_list[i]))

    ke_index  = a.lut[125]  # Kinetic Energy (KE)
    rke_index = a.lut[126]  # KE associated with radial motion
    tke_index = a.lut[127]  # KE associated with theta motion
    pke_index = a.lut[128]  # KE associated with azimuthal motion

    #We also grab some energies associated with the mean (m=0) motions
    mke_index = a.lut[129]
    mrke_index = a.lut[130]  # KE associated with mean radial motion
    mtke_index = a.lut[131]  # KE associated with mean theta motion
    mpke_index = a.lut[132]  # KE associated with mean azimuthal motion

    # ... and the energies associated with the fluctuations (i.e., 
    # energy associated with the convective amplitudes. The mean + 
    # fluctuating energies in each component should equal the total 
    # component
    fke_index = a.lut[133]
    frke_index = a.lut[134]
    ftke_index = a.lut[135]
    fpke_index = a.lut[136]

    ntimes = a.niter
    
    for j in range(ntimes):
        times.append(a.time[j])
        iters.append(a.iters[j])

        ke.append(a.vals[j,ke_index])
        rke.append(a.vals[j,rke_index])
        tke.append(a.vals[j,tke_index])
        pke.append(a.vals[j,pke_index])

        try:
            mke.append(a.vals[j, mke_index])
        except:
            mke.append(0.)

        mrke.append(a.vals[j,mrke_index])
        mtke.append(a.vals[j,mtke_index])
        mpke.append(a.vals[j,mpke_index])

        try:
            fke.append(a.vals[j, fke_index])
            frke.append(a.vals[j,frke_index])
            ftke.append(a.vals[j,ftke_index])
            fpke.append(a.vals[j,fpke_index])
        except:
            fke.append(0.)
            frke.append(0.)
            ftke.append(0.)
            fpke.append(0.)

print ('Saving file at ' + savefile + ' ...')
np.save(savefile, (times, iters, ke, rke, tke, pke, 
    mke, mrke, mtke, mpke, fke, frke, ftke, fpke))
