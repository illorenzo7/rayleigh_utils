import numpy as np
import os, sys
from diagnostic_reading import Meridional_Slice

dirname = sys.argv[1]
radatadir = dirname + '/Meridional_Slices/'

datadir = dirname + '/data/'
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)


allfiles = os.listdir(radatadir)
nallfiles = len(allfiles)
allfiles.sort()

files []
for i in range(nallfiles):
    if (allfiles[i] > cutoff_file):
        files.append(allfiles[i])

nfiles = len(files)

mer0 = Meridional_Slice(radatadir + files[0], '')
nr = mer0.nr
nt = mer0.ntheta

vp_av = np.load(datadir + 'vavg.npy')[2]

nbins = 100

# Try to determine variable range automatically for each radius
sample_vp = mer0.vals[:,:,:,mer0.lut[3],:]
min_vp = np.min(sample_vp)
max_vp = 200.

mins = -4500.
maxs = 3500.

# allow the calling routine to pass min/max values
args = sys.argv[2:]
nargs = len(args)
caller_specified_minmax = False

for ii in range(nargs):
    if (args[ii] == '-minmax'):
        caller_specified_minmax = True
        mins_maxes = args[ii + 1].split('')
        for j in range(len(mins_maxes)):
            mins_maxes[j] = np.float(mins_maxes[j])

if (caller_specified_minmax):
    min_vr = mins_maxes[0]
    maxvr = mins_maxes[1]
    mins = mins_maxes[2]
    maxs = mins_maxes[3]

# Write down the binning structure in the data directory for later
np.save(datadir + 'bin_info_vr_s.npy', (min_vr, maxvr, mins, maxs, nt, nr, nbins_vr, nbins_s))
vr_binedges = np.linspace(min_vr, maxvr, nbins_vr+1)
vr_space = (maxvr - min_vr)/nbins_vr

s_binedges = np.linspace(mins,maxs,nbins_s + 1)
s_space = (maxs - mins)/nbins_s

vr_s_dist = np.zeros((nt,nr,nbins_vr,nbins_s))

# Keep a logfile to monitor program's progress
logname = 'logfile_vr_s_dist'
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
    
    niter = mer.niter
    nphi = mer.nphi

    # Loop through each phi-value the Meridional_Slice and each
    # time in it.
    for pindex in range(nphi):
        for tindex in range(mer.niter):
            # Compute the vr' and s' values for this phi-value and time
            vr_fluc = (mer.vals[pindex,:,:,mer.lut[1],tindex] - vr_avg)/100.
            s_fluc = mer.vals[pindex,:,:,mer.lut[64],tindex] - s_avg

            # Loop through the current slice in theta-r space,
            # binning each value of vr and s accordingly.
            
            for it in range(nt):
                for ir in range(nr):
                    # compute the local vr' and s' for  this cell in the slice
                    vr_loc = vr_fluc[it,ir]
                    s_loc = s_fluc[it,ir]
                    # bin vr', putting all values outside [min_vr,maxvr] in the
                    # closest bin in range (the outermost bins)
                    if (vr_loc < min_vr + vr_space):
                        bin_index_vr = 0
                    elif (vr_loc > maxvr - vr_space):
                        bin_index_vr = nbins_vr - 1
                    else:
                        bin_index_vr = np.floor((vr_loc - min_vr)/vr_space)
                        bin_index_vr = int(bin_index_vr)

                    # bin s' according to the same process
                    if (s_loc < mins + s_space):
                        bin_index_s = 0
                    elif (s_loc > maxs - s_space):
                        bin_index_s = nbins_s - 1
                    else:
                        bin_index_s = np.floor((s_loc - mins)/s_space)
                        bin_index_s = int(bin_index_s)
                   
                    vr_s_dist[it,ir,bin_index_vr,bin_index_s] += 1.

np.save(datadir + 'vr_s_dist.npy',vr_s_dist)
