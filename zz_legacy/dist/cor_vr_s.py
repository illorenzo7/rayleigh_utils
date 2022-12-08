import numpy as np
from diagnostic_reading import *
import sys

# Get directory of run on which to do an,,sis
# and associated sub-directories
dirname = sys.argv[1]
datadir = dirname + '/data/'

# Create 'datadir' if it doesn't exist already
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

# Get basic geometry info
(rr,tt,cost,sint,rr_depth,ri,ro,d) = np.load(datadir + 'grid_info.npy')

# Get bin structure info (s units: erg/K/g, vr units: m/s)
(minvr, maxvr, mins, maxs, nt, nr, nbins_vr, nbins_s) =\
        np.load(datadir + 'bin_info.npy')
nt = int(nt); nr = int(nr); nbins_vr = int(nbins_vr); nbins_s = int(nbins_s)

# Read in area weights for each cell in the grid
(area_weights,areas,total_area,tt_weights,rr_weights) =\
    np.load(datadir + 'area_weights.npy')

# Read in the 4d vr/s distribution as fct of (r,theta)
print ('Reading in the distribution...')
vr_s_dist = np.load(datadir + 'vr_s_dist.npy')

# Normalize the distribution
print ('Normalizing the distribution...')
for it in range(nt):
    for ir in range(nr):
        norm_factor = np.sum(np.sum(vr_s_dist[it,ir,:,:]))
        vr_s_dist[it,ir,:,:] /= norm_factor


# For each bin, associate the variable with its value at bin center
dvr = (maxvr - minvr)/nbins_vr
ds = (maxs - mins)/nbins_s

vrvals = np.arange(minvr + dvr/2., maxvr + dvr/2., dvr)
svals = np.arange(mins + ds/2., maxs + ds/2., ds) 

# Compute vr values associated with upflows/downflows separately

# upflow <--> vr < 0, downflow <--> vr > 0
vrvals_df = np.arange(minvr + dvr/2., 0., dvr)
vrvals_uf = np.arange(maxvr - dvr/2., 0., -dvr)
vrvals_uf = vrvals_uf[::-1]

# Compute the probability distribution for the downflows only:
print('Computing dist(theta, r, vr, s) for downflows only...')
nvr_df = len(vrvals_df)
norm_df = np.sum(np.sum(vr_s_dist[:,:,:nvr_df,:],axis=2),axis=2)
vrs_df = vr_s_dist[:,:,:nvr_df,:]/norm_df.reshape((nt, nr, 1, 1))
# "Euthanize" vrs_df to remove nans
vrs_df = np.nan_to_num(vrs_df)

print ('...for upflows only...')
nvr_uf = len(vrvals_uf)
norm_uf = np.sum(np.sum(vr_s_dist[:,:,nvr_df:,:],axis=2),axis=2)
vrs_uf = vr_s_dist[:,:,nvr_df:,:]/norm_uf.reshape((nt, nr, 1, 1))
# Euthanize vrs_uf
vrs_uf = np.nan_to_num(vrs_uf)

# compute the correlation coefficient
print ('computing the 2d (r, theta) correlation coefficient ...')
vrvals_2d, svals_2d = np.meshgrid(vrvals, svals, indexing='ij')
vrsvals = vrvals_2d*svals_2d

cor_coef = np.zeros((nt, nr))

for ir in range(nr):
    for it in range(nt):
        dist_rt = vr_s_dist[it, ir, :, :]
        cor_vrs = np.sum(np.sum(vrsvals * dist_rt))
        vr2 = np.sum(np.sum(vrvals_2d**2 * dist_rt))
        s2 = np.sum(np.sum(svals_2d**2 * dist_rt))
        cor_coef[it, ir] = cor_vrs/np.sqrt(vr2*s2)
cor_coef = np.nan_to_num(cor_coef) # Euthanize cor_coef

# compute the correlation coefficients for upflows/downflows separately

# downflows
print ('computing correlation for downflows only...')
vrvals_2d_df, svals_2d_df = np.meshgrid(vrvals_df, svals, indexing='ij')
vrsvals_df = vrvals_2d_df*svals_2d_df

cor_coef_df = np.zeros((nt, nr))

for ir in range(nr):
    for it in range(nt):
        dist_df = vrs_df[it, ir, :, :]
        cor_df = np.sum(np.sum(vrsvals_df * dist_df))
        vr2_df = np.sum(np.sum(vrvals_2d_df**2 * dist_df))
        s2_df = np.sum(np.sum(svals_2d_df**2 * dist_df))
        cor_coef_df[it, ir] = cor_df/np.sqrt(vr2_df*s2_df)
cor_coef_df = np.nan_to_num(cor_coef_df) # Euthanize cor_coef_df

# upflows
print ('computing correlation for upflows only...')
vrvals_2d_uf, svals_2d_uf = np.meshgrid(vrvals_uf, svals, indexing='ij')
vrsvals_uf = vrvals_2d_uf*svals_2d_uf

cor_coef_uf = np.zeros((nt, nr))

for ir in range(nr):
    for it in range(nt):
        dist_uf = vrs_uf[it, ir, :, :]
        cor_uf = np.sum(np.sum(vrsvals_uf * dist_uf))
        vr2_uf = np.sum(np.sum(vrvals_2d_uf**2 * dist_uf))
        s2_uf = np.sum(np.sum(svals_2d_uf**2 * dist_uf))
        cor_coef_uf[it, ir] = cor_uf/np.sqrt(vr2_uf*s2_uf)
cor_coef_uf = np.nan_to_num(cor_coef_uf) # Euthanize cor_coef_uf

# Compute marginal correlations (functions of radius/latitude)
print('Computing marginal correlations (correlation of total flow, upflows')
print('and downflows as functions of radius and latitude separately...')

# Take into account weights of the different cells in the meridional plane
rrw_2d = rr_weights.reshape((1, nr))
ttw_2d = tt_weights.reshape((nt, 1))

cor_coef_r = np.sum(cor_coef*ttw_2d, axis=0)
cor_coef_t = np.sum(cor_coef*rrw_2d, axis=1)

ir_start = 1
ir_end = nr - 1 # ignore problematic "nan" endpoints
cor_coef_df_r = np.sum(cor_coef_df*ttw_2d, axis=0)
cor_coef_df_t = np.sum(cor_coef_df*rrw_2d, axis=1)

cor_coef_uf_r = np.sum(cor_coef_uf*ttw_2d, axis=0)
cor_coef_uf_t = np.sum(cor_coef_uf*rrw_2d, axis=1)

# Save everything
print ('Saving everything...')
np.save(datadir + 'cor_vr_s.npy', cor_coef)

np.save(datadir + 'cor_vr_s_df.npy', cor_coef_df)
np.save(datadir + 'cor_vr_s_uf.npy', cor_coef_uf)

np.save(datadir + 'cor_vr_s_r.npy', cor_coef_r)
np.save(datadir + 'cor_vr_s_t.npy', cor_coef_t)

np.save(datadir + 'cor_vr_s_df_r.npy', cor_coef_df_r)
np.save(datadir + 'cor_vr_s_df_t.npy', cor_coef_df_t)

np.save(datadir + 'cor_vr_s_uf_r.npy', cor_coef_uf_r)
np.save(datadir + 'cor_vr_s_uf_t.npy', cor_coef_uf_t)

print ('... done')
