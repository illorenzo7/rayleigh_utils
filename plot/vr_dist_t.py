import numpy as np
import matplotlib.pyplot as plt
from diagnostic_reading import *
import sys

# Get directory of run on which to do an,,sis
# and associated sub-directories
dirname = sys.argv[1]
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'

# Set defaults
saveplot = False
savefile = plotdir + 'vr_dist_t.png'
use_one_depth = False

# Read in command-line arguments (clas)
clas = sys.argv[2:]
nclas = len(clas)
lats = np.array([0., 15., 30., 45., 60., 75., 89.])
nlats = len(lats)

for ii in range(nclas):
    cla = clas[ii]
    if (cla == '-save'):
        saveplot = True
    elif (cla == '-savefile'):
        savefile = plotdir + clas[ii + 1] + '.png'
    elif (cla == '-lats'):
        lats = clas[ii + 1].split()
        nlats = len(lats)
        for jj in range(nlats):
            lats[jj] = float(lats[jj])
        lats = np.array(lats)
    elif (cla == '-depth'):
        use_one_depth = True
        depth = float(clas[ii + 1])

# Create 'plotdir' if it doesn't exist already
if (not os.path.isdir(plotdir)):
        os.makedirs(plotdir)

# Get basic geometry info
(rr,tt,cost,sint,rr_depth,ri,ro,d) = np.load(datadir + 'grid_info.npy')
latvals = (np.pi/2. - tt)*180./np.pi

# Get dist(t, vr)
if (not use_one_depth):
    vr_dist_t = np.load(datadir + 'vr_dist_t.npy')
else:
    ir = np.argmin(np.abs(rr_depth - depth))
    vr_dist = np.load(datadir + 'vr_dist.npy')
    vr_dist_t = vr_dist[:,ir,:]

# Get bin structure info (s units: erg/K/g, vr units: m/s)
(minvr, maxvr, mins, maxs, nt, nr, nbins_vr, nbins_s) =\
        np.load(datadir + 'bin_info.npy')
nt = int(nt); nr = int(nr); nbins_vr = int(nbins_vr); nbins_s = int(nbins_s)

# For each bin, associate the variable with its value at bin center
dvr = (maxvr - minvr)/nbins_vr
vrvals = np.arange(minvr + dvr/2., maxvr + dvr/2., dvr) 


# Plot the distribution at each desired latitude
for jj in range(nlats):
    it = np.argmin(np.abs(lats[jj] - latvals))
    actual_lat = latvals[it]
    dist = vr_dist_t[it, :]
    # Must take into account that 'dist' refers to the probability
    # that s' falls into a bin of width ds > 1 erg/K/g
    dist /= dvr

    plt.plot(vrvals, dist, label=r'$\lambda = %.1f^\circ$' %actual_lat)

plt.xlabel("vr' = vr - <vr> (m/s) ")
plt.ylabel('pdf ' + r'$((m/s)^{-1})$')
if (use_one_depth):
    plt.title('depth = %.1f%%' %(depth*100.))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.ylim((0,1))
plt.xlim((-300,200))
plt.legend()

if (not saveplot):
    plt.show()
else:
    plt.savefig(savefile, dpi=300)
    plt.close()
