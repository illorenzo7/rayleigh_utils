import numpy as np
import sys, os
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import colors
import matplotlib.pyplot as plt
from diagnostic_reading import Meridional_Slice
from common import strip_dirname, get_widest_range_file, get_iters_from_file

# Get the name of the run directory
dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Get grid info (if it's not computed already using grid_info.py, this will fail)
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr = len(rr)
nt = len(tt)

# First, get the correlation from the distribution at low latitudes

# Read in distribution
dist_file = get_widest_range_file(datadir, 'vr_vp_dist_full_lowlat')
print ('About to read in distribution from ' + datadir + dist_file + ' ...')
dist = np.load(datadir + dist_file)

# Get an appropriate name to save the data, fname tagged with "cor_vrvp"
iter1_dist, iter2_dist = get_iters_from_file(dist_file)
savename_dist = dirname_stripped + '_cor_vr_vp_dist_' +\
    str(iter1_dist).zfill(8) + '_' + str(iter2_dist).zfill(8) + '.npy'
    
# Get vavg file
vavg_file = get_widest_range_file(datadir, 'vavg')
vr_av, vt_av, vp_av = np.load(datadir + vavg_file)

# Average vp over low latitudes for consistency
tt_lat = tt*180/np.pi - 90
it1, it2 = np.argmin(np.abs(tt_lat - 15)), np.argmin(np.abs(tt_lat + 15))
vp_av_lowlat = np.mean(vp_av[it1:it2, :], axis=0)

# Get the bin structure for the distribution
vr_bincenters, vp_bincenters = np.load(datadir + 'vrvp_bincenters.npy')

vr_vp_fromdist = np.zeros(nr)
vr2_fromdist = np.zeros(nr)
vp2_fromdist = np.zeros(nr)

vr_vp_fromdist_df = np.zeros(nr)
vr2_fromdist_df = np.zeros(nr)
vp2_fromdist_df = np.zeros(nr)

vr_vp_fromdist_uf = np.zeros(nr)
vr2_fromdist_uf = np.zeros(nr)
vp2_fromdist_uf = np.zeros(nr)

for ir in range(nr):
#for ir in [1]:
    vp_av_loc = vp_av_lowlat[ir]/100.
    dist_loc = dist[ir]
    prob_loc = dist_loc/np.sum(dist_loc)
    vrvals, vpvals = np.meshgrid(vr_bincenters[ir], vp_bincenters[ir] - vp_av_loc, indexing='ij')
    
    # Get the probability distribution for just the downflows
    where_downflow = np.where(vrvals < 0)
    dist_df = np.zeros_like(dist_loc)
    dist_df[where_downflow] = dist_loc[where_downflow]
    prob_df = dist_df/np.sum(dist_df)
    
    # and just the upflows
    where_upflow = np.where(vrvals > 0)
    dist_uf = np.zeros_like(dist_loc)
    dist_uf[where_upflow] = dist_loc[where_upflow]
    prob_uf = dist_uf/np.sum(dist_uf)
    
#    plt.pcolormesh(vrvals, vpvals, dist_loc)
    # Average of correlations for total quantities
    vr_vp_fromdist[ir] = np.sum(vrvals*vpvals*prob_loc)
    vr2_fromdist[ir] = np.sum(vrvals**2*prob_loc)
    vp2_fromdist[ir] = np.sum(vpvals**2*prob_loc)
    
    # ...average for downflow-only quantitites
    vr_vp_fromdist_df[ir] = np.sum(vrvals*vpvals*prob_df)
    vr2_fromdist_df[ir] = np.sum(vrvals**2*prob_df)
    vp2_fromdist_df[ir] = np.sum(vpvals**2*prob_df)
    
    # ...average for downflow-only quantitites
    vr_vp_fromdist_uf[ir] = np.sum(vrvals*vpvals*prob_uf)
    vr2_fromdist_uf[ir] = np.sum(vrvals**2*prob_uf)
    vp2_fromdist_uf[ir] = np.sum(vpvals**2*prob_uf)

cor_vr_vp_fromdist = vr_vp_fromdist/np.sqrt(vr2_fromdist*vp2_fromdist)
cor_vr_vp_fromdist_df = vr_vp_fromdist_df/np.sqrt(vr2_fromdist_df*vp2_fromdist_df)
cor_vr_vp_fromdist_uf = vr_vp_fromdist_uf/np.sqrt(vr2_fromdist_uf*vp2_fromdist_uf)

# Save the correlation calculated from the distribution
print('Saving cor(vr,vp) from the distribution at '  + datadir + \
      savename_dist + ' ...')
np.save(datadir + savename_dist, (cor_vr_vp_fromdist, cor_vr_vp_fromdist_df,\
                                  cor_vr_vp_fromdist_uf))

# Now get the correlation from the actual flows (computed in the Reynolds
# stress)

# Get the Reynolds stress (must already be computed)

# Read in the Reynolds stress
rs_file = get_widest_range_file(datadir, 'rs')
print ('About to read in Reynolds stress from ' + datadir + dist_file + ' ...')
vr2_p,vt2_p,vp2_p,vrvp_p,vrvt_p,vtvp_p,\
    vr2_m,vt2_m, vp2_m, vrvp_m, vrvt_m, vtvp_m, fplus, fminus =\
    np.load(datadir + rs_file)

# Get an appropriate name to save the data, fname tagged with "cor_vrvp"
iter1_rs, iter2_rs = get_iters_from_file(rs_file)
savename_flow = dirname_stripped + '_cor_vr_vp_flow_' +\
    str(iter1_rs).zfill(8) + '_' + str(iter2_rs).zfill(8) + '.npy'


vr2_t = vr2_p + vr2_m
vp2_t = vp2_p + vp2_m
vrvp_t = vrvp_p + vrvp_m

vr2_lowlat = np.mean(vr2_t[it1:it2, :], axis=0)
vp2_lowlat = np.mean(vp2_t[it1:it2, :], axis=0)
vrvp_lowlat = np.mean(vrvp_t[it1:it2, :], axis=0)

vr2_lowlat_df = np.mean(vr2_m[it1:it2, :], axis=0)
vp2_lowlat_df = np.mean(vp2_m[it1:it2, :], axis=0)
vrvp_lowlat_df = np.mean(vrvp_m[it1:it2, :], axis=0)

vr2_lowlat_df = np.mean(vr2_m[it1:it2, :], axis=0)
vp2_lowlat_df = np.mean(vp2_m[it1:it2, :], axis=0)
vrvp_lowlat_df = np.mean(vrvp_m[it1:it2, :], axis=0)

vr2_lowlat_uf = np.mean(vr2_p[it1:it2, :], axis=0)
vp2_lowlat_uf = np.mean(vp2_p[it1:it2, :], axis=0)
vrvp_lowlat_uf = np.mean(vrvp_p[it1:it2, :], axis=0)

cor_vr_vp_fromflow = vrvp_lowlat/np.sqrt(vr2_lowlat*vp2_lowlat)
cor_vr_vp_fromflow_df = vrvp_lowlat_df/np.sqrt(vr2_lowlat_df*vp2_lowlat_df)
cor_vr_vp_fromflow_uf = vrvp_lowlat_uf/np.sqrt(vr2_lowlat_uf*vp2_lowlat_uf)

print('Saving cor(vr,vp) from the distribution at '  + datadir + \
      savename_flow + ' ...')

np.save(datadir + savename_flow, (cor_vr_vp_fromflow, cor_vr_vp_fromflow_df,\
                                  cor_vr_vp_fromflow_uf))

# Now plot the correlations
plotname_dist = dirname_stripped + '_cor_vr_vp_fromdist_' +\
    str(iter1_dist).zfill(8) + '_' + str(iter2_dist).zfill(8) + '.png'
    
plotname_flow = dirname_stripped + '_cor_vr_vp_fromflow_' +\
    str(iter1_rs).zfill(8) + '_' + str(iter2_rs).zfill(8) + '.png'  
    
plt.plot(rr/ro, cor_vr_vp_fromflow, label='tot')
plt.plot(rr/ro, cor_vr_vp_fromflow_df, label='df')
plt.plot(rr/ro, cor_vr_vp_fromflow_uf, label='uf')
plt.plot(rr/ro, np.zeros_like(rr), 'k')
plt.legend()
plt.xlabel(r'$r/r_o$', fontsize=14)
plt.ylabel(r'$cor(v_r, v_\phi)$', fontsize=14)
plt.xlim(rr[-2]/ro, rr[1]/ro)
plt.tight_layout()
# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
plt.savefig(plotdir + plotname_flow, dpi=300)
plt.close()

plt.plot(rr/ro, cor_vr_vp_fromdist, label='tot')
plt.plot(rr/ro, cor_vr_vp_fromdist_df, label='df')
plt.plot(rr/ro, cor_vr_vp_fromdist_uf, label='uf')
plt.plot(rr/ro, np.zeros_like(rr), 'k')
plt.legend()
plt.xlabel(r'$r/r_o$', fontsize=14)
plt.ylabel(r'$cor(v_r, v_\phi)$', fontsize=14)
plt.xlim(rr[-2]/ro, rr[1]/ro)
plt.tight_layout()
# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
plt.savefig(plotdir + plotname_dist, dpi=300)
plt.close()