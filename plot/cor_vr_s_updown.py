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
# Read in command-line arguments (clas)
clas = sys.argv[2:]
nclas = len(clas)
for ii in range(nclas):
    cla = clas[ii]
    if (cla == '-save'):
        saveplot = True

# Create 'plotdir' if it doesn't exist already

if (not os.path.isdir(plotdir)):
        os.makedirs(plotdir)

# Get basic geometry info
(rr,tt,cost,sint,rr_depth,ri,ro,d) = np.load(datadir + 'grid_info.npy')

# Get dist(r, vr, s)
vr_s_dist_r = np.load(datadir + 'vr_s_dist_r.npy')

# Get bin structure info (s units: erg/K/g, vr units: m/s)
(minvr, maxvr, mins, maxs, nt, nr, nbins_vr, nbins_s) =\
        np.load(datadir + 'bin_info.npy')
nt = int(nt); nr = int(nr); nbins_vr = int(nbins_vr); nbins_s = int(nbins_s)

# For each bin, associate the variable with its value at bin center
dvr = (maxvr - minvr)/nbins_vr
ds = (maxs - mins)/nbins_s

# Take vr/s bin values associated with downflow/upflow only 
# (vr < 0)/(vr > 0)
vrvals_df = np.arange(minvr + dvr/2., 0., dvr)
vrvals_uf = np.arange(maxvr - dvr/2., 0., -dvr)
vrvals_uf = vrvals_uf[::-1]

# Take all svals
svals = np.arange(mins + ds/2., maxs + ds/2., ds) 

# Re-normalize the probability distribution over downflows only:
nvr_df = len(vrvals_df)
nvr_uf = len(vrvals_uf)
rn_df = np.sum(np.sum(vr_s_dist_r[:,:nvr_df,:],axis=1),axis=1)
rn_uf = np.sum(np.sum(vr_s_dist_r[:,nvr_df:,:],axis=1),axis=1)

vrs_df = vr_s_dist_r[:,:nvr_df,:]/(rn_df.reshape((nr,1,1)))
vrs_uf = vr_s_dist_r[:,nvr_df:,:]/(rn_uf.reshape((nr,1,1)))



# compute the correlation coefficients for upflows/downflows separately

vrvals_2d_df, svals_2d_df = np.meshgrid(vrvals_df, svals, indexing='ij')
vrvals_2d_uf, svals_2d_uf = np.meshgrid(vrvals_uf, svals, indexing='ij')

vrsvals_df = vrvals_2d_df*svals_2d_df
vrsvals_uf = vrvals_2d_uf*svals_2d_uf

# downflows
cor_df = np.sum(np.sum(vrsvals_df.reshape((1, nvr_df, nbins_s))\
        *vrs_df, axis=1), axis=1)
vr2_df = np.sum(np.sum(vrvals_2d_df.reshape((1,nvr_df,nbins_s))**2\
        *vrs_df, axis=1), axis=1)

s2_df = np.sum(np.sum(svals_2d_df.reshape((1,nvr_df,nbins_s))**2\
        *vrs_df, axis=1), axis=1)

cor_coef_df = cor_df/np.sqrt(vr2_df*s2_df)

# upflows
cor_uf = np.sum(np.sum(vrsvals_uf.reshape((1, nvr_uf, nbins_s))\
        *vrs_uf, axis=1), axis=1)
vr2_uf = np.sum(np.sum(vrvals_2d_uf.reshape((1,nvr_uf,nbins_s))**2\
        *vrs_uf, axis=1), axis=1)

s2_uf = np.sum(np.sum(svals_2d_uf.reshape((1,nvr_uf,nbins_s))**2\
        *vrs_uf, axis=1), axis=1)

cor_coef_uf = cor_uf/np.sqrt(vr2_uf*s2_uf)


plt.plot(rr_depth[1:],cor_coef_df[1:],color='r',label='downflow')
plt.plot(rr_depth[1:],cor_coef_uf[1:],color='b',label='upflow')

plt.xlabel('depth = (ro - r)/(ro - ri)')
plt.ylabel("<vr's'>/sqrt(<(vr')^2> <(s')^2>)")
plt.plot(.05*np.ones(100),np.linspace(0,1,100),'k--')
plt.ylim((0,1))
plt.xlim((0,1))
plt.legend()

if (not saveplot):
    plt.show()
else:
    plt.savefig(plotdir + 'cor_vr_s_updown.png', dpi=300)
    plt.close()
