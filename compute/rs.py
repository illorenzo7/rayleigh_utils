# Author: Loren Matilsky
# Created: 04/26/2018
# Comments added: 05/07/2018
# This script computes all six components of the Reynolds stress tensor due 
# to the upflows and downflows separately for a Rayleigh run (stored in dir-
# name.) The script uses the run's Meridional Slices to get instantaneous 
# snapshots of the flow in the meridional plane and separate the various 
# correlations into their contributions due to upflows and downflows. Filling
# factors for the upflows and downflows (relative number of times an upflow
# or downflow is recorded at a particular point in the meridional plane) are
# also computed. The output is stored as
# (vr2_p, vt2_p, vp2_p, vrvp_p, vrvt_p, vtvp_p,
#  vr2_m, vt2_m, vp2_m, vrvp_m, vrvt_m, vtvp_m, fplus, fminus)
# in the file [dirname]_rs_[first iter]_[last iter].npy
# Units are in (cm/s)^2
# Quantities associated with v_phi have v_phi's deviation from the azimuthal
# mean (i.e., the differential rotation has been subtracted out)

# Import relevant modules
import numpy as np
import os, sys
from diagnostic_reading import Meridional_Slice
from common import get_file_lists, get_desired_range, strip_dirname,\
    get_widest_range_file

# Function to split an array according to locations where it's positive or
    # negative
def split_pm(vr):
    index_minus = np.where(vr < 0)
    index_plus = np.where(vr > 0)
    return(index_plus, index_minus)

# Multiply two arrays while applying the mask indicated by "index_where"; 
    # values not indicated by the mask will be set to zero
def prod(v1, v2, index_where):
    v1v2 = v1*v2   
    v1v2_where = np.zeros_like(v1v2)
    v1v2_where[index_where] = v1v2[index_where]
    return v1v2_where

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
data_type = 'Meridional_Slices'
radatadir = dirname + '/' + data_type + '/'

# Get full list of files in radatadir
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Get grid info (if it's not computed already using grid_info.py, this will fail)
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr = len(rr)
nt = len(tt)

# Read in CLAs (if any) to change default file range to use (by default, this
# will be the last 100 meridional slices)
args = sys.argv[2:]
nargs = len(args)
if (nargs == 0):
    index_first, index_last = nfiles - 100, nfiles
else:
    index_first, index_last = get_desired_range(int_file_list, args)

print ('Making distribution for files in the range %s to %s' 
       %(file_list[index_first], file_list[index_last - 1]))
   
# Set the savename by the directory, what we are saving, and first and last
# iteration files for the average
savename = dirname_stripped + '_rs_' + file_list[index_first] + '_' +\
    file_list[index_last - 1] + '.npy'
savefile = datadir + savename  

# Get vavg file and data
vavg_file = get_widest_range_file(datadir, 'vavg')
print ('getting vavg from %s...' %(datadir + vavg_file))
vr_av, vt_av, vp_av = np.load(datadir + vavg_file)

# Filling factors for upflows and downflows
fplus = np.zeros_like(vp_av)
fminus = np.zeros_like(vp_av)

vr2_p = np.zeros_like(vp_av)
vr2_m = np.zeros_like(vp_av)

vt2_p = np.zeros_like(vp_av)
vt2_m = np.zeros_like(vp_av)

vp2_p = np.zeros_like(vp_av)
vp2_m = np.zeros_like(vp_av)


vrvp_p = np.zeros_like(vp_av)
vrvp_m = np.zeros_like(vp_av)

vrvt_p = np.zeros_like(vp_av)
vrvt_m = np.zeros_like(vp_av)

vtvp_p = np.zeros_like(vp_av)
vtvp_m = np.zeros_like(vp_av)

count = 0
for ii in range(index_first, index_last):
    mer = Meridional_Slice(radatadir + file_list[ii], '')
    print ('considering Meridional_Slices/%s for the Reynolds stress ...'\
           %file_list[ii])
    niter = mer.niter
    nphi = mer.nphi

    # Loop through each phi-value the Meridional_Slice and each
    # time in it.
    for tindex in range(mer.niter):
        for pindex in range(nphi):
            vr_merslice = mer.vals[pindex, :, :, mer.lut[1], tindex]
            vt_merslice = mer.vals[pindex, :, :, mer.lut[2], tindex]
            vp_merslice_fluc = mer.vals[pindex, :, :, mer.lut[3], tindex] - vp_av

            # Split things according to upflows and downflows
            index_p, index_m = split_pm(vr_merslice)

            fplus[index_p] += 1.
            fminus[index_m] += 1.

            vr2_p += prod(vr_merslice, vr_merslice, index_p)
            vr2_m += prod(vr_merslice, vr_merslice, index_m)

            vt2_p += prod(vt_merslice, vt_merslice,index_p)
            vt2_m += prod(vt_merslice, vt_merslice,index_m)

            vp2_p += prod(vp_merslice_fluc, vp_merslice_fluc,index_p)
            vp2_m += prod(vp_merslice_fluc, vp_merslice_fluc,index_m)

            vrvp_p += prod(vr_merslice, vp_merslice_fluc,index_p)
            vrvp_m += prod(vr_merslice, vp_merslice_fluc,index_m)

            vrvt_p += prod(vr_merslice, vt_merslice,index_p)
            vrvt_m += prod(vr_merslice, vt_merslice,index_m)

            vtvp_p += prod(vt_merslice, vp_merslice_fluc,index_p)
            vtvp_m += prod(vt_merslice, vp_merslice_fluc,index_m)
            count += 1

vr2_m /= count
vr2_p /= count
vt2_m /= count
vt2_p /= count
vp2_m /= count
vp2_p /= count

vrvp_m /= count
vrvp_p /= count
vrvt_p /= count
vrvt_m /= count
vtvp_m /= count
vtvp_p /= count

fplus /= count
fminus /= count

print ('Saving the Reynolds stresses at ' + datadir + savename + ' ...')
np.save(datadir + savename,(vr2_p,vt2_p,vp2_p,vrvp_p,vrvt_p,vtvp_p,\
    vr2_m,vt2_m,vp2_m, vrvp_m, vrvt_m, vtvp_m, fplus, fminus))