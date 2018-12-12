import sys
import numpy as np
from common import strip_dirname, get_widest_range_file, get_iters_from_file
# Get the name of the run directory
dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist
datadir = dirname + '/data/'

# Get location of thermal BL 
ir_tbl = np.load(datadir + dirname_stripped + '_ir_tbl.npy')

# Get the downflow distribution from the pre-computed "split" files
################################################################
# 1, lowlat, or 0 to 15 degrees latitude North and South
dist_file = get_widest_range_file(datadir, 'vr_s_dist_df_1_lowlat')
iter1, iter2 = get_iters_from_file(dist_file)
print ('About to read in distribution from ' + datadir + dist_file + ' ...')
dist = np.load(datadir + dist_file)
dist_tbl = dist[ir_tbl]
dist_tbl[-1, :] = dist_tbl[-2, :] # "Fix the fuck-up I made on the last bin
savename = dirname_stripped + '_vr_s_dist_df_1_lowlat_tbl_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.npy'
savefile = datadir + savename
print ('Saving ' + savefile + ' ...')
np.save(savefile, dist_tbl)

# 2, lowmidlat, or 15 to 30 degrees latitude North and South
dist_file = get_widest_range_file(datadir, 'vr_s_dist_df_2_lowmidlat')
iter1, iter2 = get_iters_from_file(dist_file)
print ('About to read in distribution from ' + datadir + dist_file + ' ...')
dist = np.load(datadir + dist_file)
dist_tbl = dist[ir_tbl]
dist_tbl[-1, :] = dist_tbl[-2, :] # "Fix the fuck-up I made on the last bin
savename = dirname_stripped + '_vr_s_dist_df_2_lowmidlat_tbl_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.npy'
savefile = datadir + savename
print ('Saving ' + savefile + ' ...')
np.save(savefile, dist_tbl)

# 3, midlat, or 30 to 45 degrees latitude North and South
dist_file = get_widest_range_file(datadir, 'vr_s_dist_df_3_midlat')
iter1, iter2 = get_iters_from_file(dist_file)
print ('About to read in distribution from ' + datadir + dist_file + ' ...')
dist = np.load(datadir + dist_file)
dist_tbl = dist[ir_tbl]
dist_tbl[-1, :] = dist_tbl[-2, :] # "Fix the fuck-up I made on the last bin
savename = dirname_stripped + '_vr_s_dist_df_3_midlat_tbl_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.npy'
savefile = datadir + savename
print ('Saving ' + savefile + ' ...')
np.save(savefile, dist_tbl)

# 4, midhighlat, or 45 to 60 degrees latitude North and South
dist_file = get_widest_range_file(datadir, 'vr_s_dist_df_4_midhighlat')
iter1, iter2 = get_iters_from_file(dist_file)
print ('About to read in distribution from ' + datadir + dist_file + ' ...')
dist = np.load(datadir + dist_file)
dist_tbl = dist[ir_tbl]
dist_tbl[-1, :] = dist_tbl[-2, :] # "Fix the fuck-up I made on the last bin
savename = dirname_stripped + '_vr_s_dist_df_4_midhighlat_tbl_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.npy'
savefile = datadir + savename
print ('Saving ' + savefile + ' ...')
np.save(savefile, dist_tbl)

# 5, highlat, or 60 to 75 degrees latitude North and South
dist_file = get_widest_range_file(datadir, 'vr_s_dist_df_5_highlat')
iter1, iter2 = get_iters_from_file(dist_file)
print ('About to read in distribution from ' + datadir + dist_file + ' ...')
dist = np.load(datadir + dist_file)
dist_tbl = dist[ir_tbl]
dist_tbl[-1, :] = dist_tbl[-2, :] # "Fix the fuck-up I made on the last bin
savename = dirname_stripped + '_vr_s_dist_df_5_highlat_tbl_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.npy'
savefile = datadir + savename
print ('Saving ' + savefile + ' ...')
np.save(savefile, dist_tbl)

# 6, superhighlat, or 75 to 90 degrees latitude North and South
dist_file = get_widest_range_file(datadir, 'vr_s_dist_df_6_superhighlat')
iter1, iter2 = get_iters_from_file(dist_file)
print ('About to read in distribution from ' + datadir + dist_file + ' ...')
dist = np.load(datadir + dist_file)
dist_tbl = dist[ir_tbl]
dist_tbl[-1, :] = dist_tbl[-2, :] # "Fix the fuck-up I made on the last bin
savename = dirname_stripped + '_vr_s_dist_df_6_superhighlat_tbl_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.npy'
savefile = datadir + savename
print ('Saving ' + savefile + ' ...')
np.save(savefile, dist_tbl)