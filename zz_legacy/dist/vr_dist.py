import numpy as np
import os, sys
from diagnostic_reading import Meridional_Slice

dirname = sys.argv[1]
radatadir = dirname + '/Meridional_Slices/'

datadir = dirname + '/data/'
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)


files = os.listdir(radatadir)
nfiles = len(files)
files.sort()

#times,dtimes,weights = np.load(datadir + 'mer_times.npy')

mer0 = Meridional_Slice(radatadir + files[0], '')
nr = mer0.nr
nt = mer0.ntheta

vavg = np.load(datadir + 'vavg.npy')
vr_av = vavg[0]

# make 100 bins for vr in the range [-300 m/s, 200 m/s]
nbins = 100
binedges = np.linspace(-300., 200., nbins+1)
minvr = np.min(binedges)
maxvr = np.max(binedges)

counts = np.zeros((nt,nr,nbins))

for ii in range(nfiles):
#for ii in range(10): # for debugging purposes
    mer = Meridional_Slice(radatadir + files[ii], '')
    print ('adding mer slice %s to the average...' %files[ii])
    niter = mer.niter
    nphi = mer.nphi

    for tindex in range(mer.niter):
        for pindex in range(nphi):
            vr_fluc = (mer.vals[pindex,:,:,mer.lut[1],tindex] - vr_av)/100.
            # points outside the total bin range add to the closest bin
            # inside the range:
            indicator_low = np.zeros_like(vr_av)
            indicator_low[np.where(vr_fluc < minvr)] = 1.
            counts[:,:,0] += indicator_low

            indicator_high = np.zeros_like(vr_av)
            indicator_high[np.where(vr_fluc >= maxvr)] = 1.
            counts[:,:,nbins-1] += indicator_high
            
            # now loop through all the bins in range:
            for ii in range(nbins):
                leftedge = binedges[ii]
                rightedge = binedges[ii+1]
                indicator = np.zeros_like(vr_av)
                indicator[np.where((vr_fluc >= leftedge)*(vr_fluc < rightedge))] = 1.
                counts[:,:,ii] += indicator


np.save(datadir + 'counts.npy',counts)


