import numpy as np
import os, sys

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]

# Directory with data
datadir = dirname + '/data/'

# Create 'datadir' if it doesn't exist already
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

# Read in area weights for each cell in the grid
(area_weights,areas,total_area,tt_weights,rr_weights)\
        = np.load(datadir + 'area_weights.npy')

# Read in the 4d vr/s distribution as fct of (r,theta)
print ('Reading in the distribution...')
vr_s_dist = np.load(datadir + 'vr_s_dist.npy')

# Shape of distribution
nt, nr, vr_nbins, s_nbins = np.shape(vr_s_dist)

# Normalize the distribution
print ('Normalizing the distribution...')
for it in range(nt):
    for ir in range(nr):
        vr_s_dist[it,ir,:,:] /= \
                (np.sum(np.sum(vr_s_dist[it,ir,:,:], axis=0), axis=0))

# Weights in 3d and 4d to marginalize the 4d distribution
rrw_3d = rr_weights.reshape((1,nr,1))
ttw_3d = tt_weights.reshape((nt,1,1))
rrw_4d = rr_weights.reshape((1,nr,1,1))
ttw_4d = tt_weights.reshape((nt,1,1,1))


print ('Marginalizing dist(theta, r, vr, s) over theta and r...')
# Marginalize dist over theta --> dist(r, vr, s)
vr_s_dist_r = np.sum(ttw_4d*vr_s_dist,axis=0)
# ...over r --> dist(theta, vr, s)
vr_s_dist_t = np.sum(rrw_4d*vr_s_dist,axis=1)

print ('Marginalizing dist(theta, r, vr, s) over s and vr')
print ('         --> vr_dist and s_dist...')
# ...over s --> vr_dist := dist(theta, r, vr)
vr_dist = np.sum(vr_s_dist,axis=3)
# ...over vr --> s_dist := dist(theta, r, s)
s_dist = np.sum(vr_s_dist, axis=2)

print ('Marginalizing vr_dist over theta and r...')
# Marginalize vr_dist over theta --> dist(r, vr)
vr_dist_r = np.sum(vr_dist*ttw_3d, axis=0)
# ...vr_dist over r --> dist(theta, vr)
vr_dist_t = np.sum(vr_dist*rrw_3d, axis=1)

print ('Marginalizing s_dist over theta and r...')
# Marginalize s_dist over theta --> dist(r, s)
s_dist_r = np.sum(s_dist*ttw_3d, axis=0)
# ...over r --> dist(theta, s)
s_dist_t = np.sum(s_dist*rrw_3d, axis=1)

# Save everything
print ('Saving everything...')
np.save(datadir + 'vr_s_dist_r.npy',vr_s_dist_r)
np.save(datadir + 'vr_s_dist_t.npy',vr_s_dist_t)

np.save(datadir + 'vr_dist.npy',vr_dist)
np.save(datadir + 's_dist.npy',s_dist)

np.save(datadir + 'vr_dist_r', vr_dist_r)
np.save(datadir + 'vr_dist_t', vr_dist_t)

np.save(datadir + 's_dist_r', s_dist_r)
np.save(datadir + 's_dist_t', s_dist_t)
print ('...done')
