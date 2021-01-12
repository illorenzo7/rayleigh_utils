# Author: Loren Matilsky
# Created: 02/25/2018

# This script computes the 2-D distribution for v_r and v_phi (radial and 
# azimuthal velocities at each point in the N_theta x N_r meridional plane. 

# Takes the main directory of a Rayleigh run as input. 

# The script makes use of the Meridional_Slice object in rayleigh, which 
# has the following attributes:

#    self.niter                                    : number of time steps
#    self.nq                                       : number of diagnostic quantities output
#    self.nr                                       : number of radial points
#    self.ntheta                                   : number of theta points
#    self.nphi                                     : number of phi points sampled
#    self.qv[0:nq-1]                               : quantity codes for the diagnostics output
#    self.radius[0:nr-1]                           : radial grid
#    self.costheta[0:ntheta-1]                     : cos(theta grid)
#    self.sintheta[0:ntheta-1]                     : sin(theta grid)
#    self.phi[0:nphi-1]                            : phi values (radians)
#    self.phi_indices[0:nphi-1]                    : phi indices (from 1 to nphi)
#    self.vals[0:nphi-1,0:ntheta-1,0:nr-1,0:nq-1,0:niter-1] : The meridional slices 
#    self.iters[0:niter-1]                         : The time step numbers stored in this output file
#    self.time[0:niter-1]                          : The simulation time corresponding to each time step
#    self.version                                  : The version code for this particular output (internal use)
#    self.lut                                      : Lookup table for the different diagnostics output

import numpy as np
import os, sys
from diagnostic_reading import Meridional_Slice

# Read in the run directory
dirname = sys.argv[1]
radatadir = dirname + '/Meridional_Slices/'

datadir = dirname + '/data/'
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

# Read in list of meridional slice files
files = os.listdir(radatadir)
nfiles = len(files)
files.sort()

# Get basic grid structure from first meridional slice
mer0 = Meridional_Slice(radatadir + files[0], '')
nr = mer0.nr
nt = mer0.ntheta

# make many bins for vr and vp
nbins_vr = 100
nbins_vp = 100

# Array to hold the distribution
vr_vp_dist = np.zeros((nt, nr, nbins_vr, nbins_vp))

# Read in the bin limits (function of radius; multiply by 1.05 to be safe)
# Note the limits are in m/s
vr_maxes, vr_mins, vp_maxes, vp_mins = \
        np.load(datadir + 'vr_vp_dist_limits.npy')
vr_maxes *= 1.05
vp_maxes *= 1.05
vr_mins *= 1.05
vp_mins *= 1.05


# Compute binning structure for each radius
vr_binedges = np.zeros((nr, nbins_vr + 1))
vr_bincenters = np.zeros((nr, nbins_vr))

vp_binedges = np.zeros((nr, nbins_vr + 1))
vp_bincenters = np.zeros((nr, nbins_vr))

vr_space = np.zeros(nr)
vp_space = np.zeros(nr)

def vrvp_sort(vr_mins, vr_maxes, vp_
for ir in range(nr):
    minvr = vr_mins[ir]
    maxvr = vr_maxes[ir]
    minvp = vp_mins[ir]
    maxvp = vp_maxes[ir]

    vr_binedges[ir,:] = np.linspace(minvr, maxvr, nbins_vr+1)
    vr_bincenters = 0.5*(vr_binedges[:-1] + vr_binedges[1:])
    vr_space[ir] = (maxvr - minvr)/nbins_vr

    vp_binedges[ir,:] = np.linspace(minvp, maxvp, nbins_vp+1)
    vp_bincenters = 0.5*(vp_binedges[:-1] + vp_binedges[1:])
    vp_space[ir] = (maxvp - minvp)/nbins_vp

# Save binspacing and binedges to use later in plots
np.save(datadir + 'vrvp_binedges.npy', (vr_binedges, vp_binedges))
np.save(datadir + 'vrvp_bincenters.npy', (vr_bincenters, vp_bincenters))
np.save(datadir + 'vrvp_binspacing.npy', (vr_space, vp_space))

# Keep a logfile to monitor program's progress
logname = 'logfile_vr_vp_dist'
logfile = open(datadir + logname,'w')
opened_only_once = True

for ii in range(nfiles):
#for ii in range(10): # for debugging purposes
    mer = Meridional_Slice(radatadir + files[ii], '')

    # Log the progress we are making!
    if (not opened_only_once):
        logfile = open(datadir + logname, 'a')
    else:
        opened_only_once = False
    logfile.write('adding %s/Meridional_Slices/%s to the distribution...\n' %(dirname.split('/')[-1],files[ii]))
    logfile.close()
    
    # Print the message for good measure
    print('adding %s/Meridional_Slices/%s to the distribution...\n' %(dirname.split('/')[-1],files[ii]))
    
    niter = mer.niter
    nphi = mer.nphi

    # Loop through each phi-value the Meridional_Slice and each
    # time in it.
    for pindex in range(nphi):
        for tindex in range(mer.niter):
            # Compute the vr and vp values for this phi-value and time
            # convert data from (cm/s) --> (m/s)
            vr_mer = mer.vals[pindex,:,:,mer.lut[1],tindex]/100.
            vp_mer = mer.vals[pindex,:,:,mer.lut[3],tindex]/100.

            # Loop through the current slice in theta-r space,
            # binning each value of vr and s accordingly.
            
            for ir in range(nr):

                # Get max bin values and spacing
                minvr_loc = vr_mins[ir]
                maxvr_loc = vr_maxes[ir]
                minvp_loc = vp_mins[ir]
                maxvp_loc = vp_maxes[ir]

                vr_space_loc = vr_space[ir]
                vp_space_loc = vp_space[ir]

                for it in range(nt):
                    # compute the local vr/vp for  this cell in the slice
                    vr_loc = vr_mer[it,ir]
                    vp_loc = vp_mer[it,ir]
                    # bin vr, putting all values outside [minvr,maxvr] 
                    # in the
                    # closest bin in range (the outermost bins)
                    if (vr_loc < minvr_loc + vr_space_loc):
                        bin_index_vr = 0
                    elif (vr_loc > maxvr_loc - vr_space_loc):
                        bin_index_vr = nbins_vr - 1
                    else:
                        bin_index_vr = np.floor((vr_loc - minvr_loc)\
                                /vr_space_loc)
                        bin_index_vr = int(bin_index_vr)

                    # bin vp according to the same process
                    if (vp_loc < minvp_loc + vp_space_loc):
                        bin_index_vp = 0
                    elif (vp_loc > maxvp_loc - vp_space_loc):
                        bin_index_vp = nbins_vp - 1
                    else:
                        bin_index_vp = np.floor((vp_loc - minvp_loc)\
                                /vp_space_loc)
                        bin_index_vp = int(bin_index_vp)
                   
                    vr_vp_dist[it,ir,bin_index_vr,bin_index_vp] += 1.

np.save(datadir + 'vr_vp_dist.npy',vr_vp_dist)
