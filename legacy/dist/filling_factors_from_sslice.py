import numpy as np
import matplotlib.pyplot as plt
from diagnostic_reading import ShellSlice
import sys, os

dirname = sys.argv[1]
datadir = dirname + '/data/'

# Create 'datadir' if it doesn't exist already
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

slicedir = dirname + '/Shell_Slices/'
files = os.listdir(slicedir)
files.sort()

rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr = len(rr)
nt = len(tt)
nph = 2*nt
sint_2d = sint.reshape((1, nt))

print ('Reading ' + slicedir + files[-1] + ' ...')
a = ShellSlice(slicedir + files[-1], '')

vr = a.vals[:, :, :, a.lut[1], 0]

nrad = len(a.radius)
fplus = np.zeros(nrad)
fminus = np.zeros(nrad)

for i in range(nrad):
    vr_loc = vr[:, :, i]
    where_p = np.where(vr_loc > 0)
    where_m = np.where(vr_loc < 0)
    ind_p = np.zeros((nph, nt))
    ind_p[where_p] = 1
    ind_m = np.zeros((nph, nt))
    ind_m[where_m] = 1
    
    num_p = np.sum(ind_p*sint_2d)/np.sum(sint_2d)/nph
    num_m = np.sum(ind_m*sint_2d)/np.sum(sint_2d)/nph
    fplus[i] = num_p
    fminus[i] = num_m
    
plt.scatter(a.radius/ro, fplus, label='up')
plt.scatter(a.radius/ro, fminus, label='down')
plt.title(dirname_stripped, fontsize=18)
# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
plt.legend()
plt.xlim((ri/ro, 1))
plt.ylim((0, 1))
plt.xlabel(r'$r/r_o$')
plt.ylabel('filling factor')
