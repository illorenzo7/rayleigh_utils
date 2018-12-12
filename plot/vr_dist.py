import numpy as np
import os
import sys
from diagnostic_reading import Meridional_Slice
import matplotlib.pyplot as plt

dirname = sys.argv[1]
radatadir = dirname + '/Meridional_Slices/'

datadir = dirname + '/data/'
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

files = os.listdir(radatadir)
nfiles = len(files)
files.sort()

# parse CLAs
saveplot=False
args = sys.argv[1:]
nargs = len(args)
for ii in range(nargs):
    if (args[ii] == '-depths'):
        depths = args[ii + 1].split()
        for j in range(len(depths)):
            depths[j] = float(depths[j])
    elif (args[ii] == '-save'):
        saveplot = True
    elif (args[ii] == '-savefile'):
        savefile = args[ii + 1]

#times,dtimes,weights = np.load(datadir + 'mer_times.npy')

mer0 = Meridional_Slice(radatadir + files[0], '')
nr = mer0.nr
nt = mer0.ntheta
rr = mer0.radius
ro = rr[0]
ri = rr[nr-1]
depth = ro - ri

vr_dist = np.load(datadir + 'dist_short.npy')

area_weights,areas,total_area,tt_weights,rr_weights = np.load(datadir + 'area_weights.npy')

tt_weights_3d = tt_weights.reshape((nt,1,1))

vr_dist_r = np.mean(vr_dist*tt_weights_3d, axis=0)

# get probability in inverse m/s:
norms = np.sum(vr_dist_r,axis=1)
norms = norms.reshape(nr,1)

prob = vr_dist_r/norms #now in terms of dimensionless probability
        # must further normalize to --> probability of being found in a bin
        # of width 5. m/s (100 bins)

prob /= 5. #units of probability are now (m/s)^(-1)

#depths = [0.01,0.02,0.03,0.04,0.05] #radii whose histograms we plot in units of fractional depth from the upper boundary
ndepths = len(depths)

rinds = np.zeros_like(depths,dtype=int)

for ir in range(ndepths):
    frac_depth = depths[ir]
    true_depth = frac_depth*depth
    desired_radius = ro - true_depth

    diffs = np.abs(desired_radius - rr)
    rinds[ir] = np.argmin(diffs)

vr_range = np.linspace(-300.,200.,100)

colors = ['b','g','r','m','k','y','c','darkgreen','orange','pink','brown']

for ir in range(ndepths):
    rind = rinds[ir]
    rad = rr[rind]
    frac_depth = (ro - rad)/depth 
    plt.plot(vr_range,prob[rind,:],label=('depth = %.1f%%' %(100*frac_depth)),color=colors[ir] )


plt.legend(loc='upper left')
plt.xlabel('vr (m/s)')
plt.ylabel('pdf' +r'$((m/s)^{-1})$')

if (saveplot):
    plt.savefig(datadir + savefile + '.png',dpi=300)
else:
    plt.show()
