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
savefile = plotdir + 's_dist_r.png'

# Read in command-line arguments (clas)
clas = sys.argv[2:]
nclas = len(clas)
depths = np.array([.01, .05, .1, .2, .35, .5, .75, .99])
ndepths = len(depths)

for ii in range(nclas):
    cla = clas[ii]
    if (cla == '-save'):
        saveplot = True
    elif (cla == '-savefile'):
        savefile = plotdir + clas[ii + 1] + '.png'
    elif (cla == '-depths'):
        depths = clas[ii + 1].split()
        ndepths = len(depths)
        for jj in range(ndepths):
            depths[jj] = float(depths[jj])
        depths = np.array(depths)

# Create 'plotdir' if it doesn't exist already
if (not os.path.isdir(plotdir)):
        os.makedirs(plotdir)

# Get basic geometry info
(rr,tt,cost,sint,rr_depth,ri,ro,d) = np.load(datadir + 'grid_info.npy')

# Get dist(r, vr, s)
s_dist_r = np.load(datadir + 's_dist_r.npy')

# Get bin structure info (s units: erg/K/g, vr units: m/s)
(minvr, maxvr, mins, maxs, nt, nr, nbins_vr, nbins_s) =\
        np.load(datadir + 'bin_info.npy')
nt = int(nt); nr = int(nr); nbins_vr = int(nbins_vr); nbins_s = int(nbins_s)

# For each bin, associate the variable with its value at bin center
ds = (maxs - mins)/nbins_s
svals = np.arange(mins + ds/2., maxs + ds/2., ds) 

# Get indices corresponding to depth values
ir_vals = np.zeros(ndepths, dtype=int)
for jj in range(ndepths):
    ir_vals[jj] = np.argmin(np.abs(rr_depth - depths[jj]))

# Plot the distribution at each radius:
for jj in range(ndepths):
    ir = ir_vals[jj]
    depth = 100.*(rr_depth[ir])
    dist = s_dist_r[ir, :]
    # Must take into account that 'dist' refers to the probability
    # that s' falls into a bin of width ds > 1 erg/K/g
    dist /= ds

    plt.plot(svals, dist, label='depth = %.1f%%' %depth)

plt.xlabel("s' = s - <s> (erg/K/g) ")
plt.ylabel('pdf ' + r'$((erg/K/g)^{-1})$')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.ylim((0,1))
#plt.xlim((0,1))
plt.legend()

if (not saveplot):
    plt.show()
else:
    plt.savefig(savefile, dpi=300)
    plt.close()
