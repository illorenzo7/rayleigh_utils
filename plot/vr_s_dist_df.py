# Author: Loren Matilsky
# Created: 05/07/2018
# This script takes output from compute/vr_s_df_split_in_latitude.py and creates
# plots at each radial coordinate of the joint v_r and S' distribution. nr 
# plots each at six different latitude ranges are produced, and saved in 
# the following subdirectories of plots/dist/vr_s_df/:
# 1_lowlat (+/- 15 degrees)
# 2_lowmidlat (15 to 30 degrees North and South)
# 3_midlat (30 to 45 degrees North and South)
# 4_midhighlat (45 to 60 degrees North and South)
# 5_highlat (60 to 75 degrees North and South)
# 6_superhighlat (75 to 90 degrees North and South)

import numpy as np
import sys, os
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import colors
import matplotlib.pyplot as plt
from common import strip_dirname, get_widest_range_file

# Get the name of the run directory
dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)

# Find the relevant places to store the plots and create them if they
# don't already exist
datadir = dirname + '/data/'
base_plotdir = dirname + '/plots/dist/vr_s_df/'
plotdirs = [base_plotdir + '1_lowlat/', base_plotdir + '2_lowmidlat/',\
            base_plotdir + '3_midlat/', base_plotdir + '4_midhighlat/',\
            base_plotdir + '5_highlat/', base_plotdir + '6_superhighlat/']
for plotdir in plotdirs:
    if (not os.path.isdir(plotdir)):
        os.makedirs(plotdir)

# Get grid info (if it's not computed already using grid_info.py, this will fail)
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr = len(rr)
nt = len(tt)

# Get the bin structure for the distribution
vr_bincenters, s_bincenters = np.load(datadir + 'vrs_df_bincenters.npy')
nbins_vr, nbins_s = len(vr_bincenters[0]), len(s_bincenters[0])

# Run a plotting loop (1 plot for each radius) at each of the six different
# latitude ranges

# Read in distribution
dist_file = get_widest_range_file(datadir, 'vr_s_dist_df_1_lowlat')
print ('About to read in distribution from ' + datadir + dist_file + ' ...')
dist = np.load(datadir + dist_file)

# Set linewidth
lw = 1

cbar_maxes = np.zeros(nr)
cbar_mins = np.zeros(nr)

#if (True): # Debugging purposes
for i in range(nr):
#    i=1 # Debugging purposes
    dist_depth = dist[i]
    cbar_mins[i] = .1
    cbar_maxes[i] = dist_depth.max()
    depth = rr_depth[i]
    # Bin values for this depth
    vrvals, svals = np.meshgrid(vr_bincenters[i], s_bincenters[i])    
    plt.pcolormesh(vrvals, svals, np.transpose(dist_depth), cmap='Greens',\
                   norm=colors.LogNorm(vmin=.1, vmax=dist_depth.max()))
    plt.plot(vr_bincenters[i], np.zeros(nbins_vr), 'k', linewidth=lw)
        # Mark line S' = 0
    plt.xlabel(r'$v_r\ (\rm{m/s})$', fontsize=14)
    plt.ylabel(r'$(S - \langle S\rangle_s)\ (\rm{erg/K/g})$',fontsize=14)
#    plt.xlim((np.min(vrvals), 0))
    plt.title(dirname_stripped + r'$\ v_r/S^\prime$' +\
              ' dist (df, lowlat), depth = %.2f%%' %(depth*100), fontsize=14)
    # Get ticks everywhere
    plt.minorticks_on()
    plt.tick_params(top='on', right='on', direction='in', which='both')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('counts', rotation=270)
    plt.tight_layout()
    savename = dirname_stripped + '_vr_s_df_depth' +\
        str(i).zfill(3) + '.png'
    print('saving plot at ' + plotdirs[0] + savename + ' ...')
    plt.savefig(plotdirs[0] + savename, dpi=300)
    plt.close()
np.save(datadir + 'vrs_df_cbar_limits_00_to_15.npy', (cbar_maxes, cbar_mins))
  
# Read in distribution
dist_file = get_widest_range_file(datadir, 'vr_s_dist_df_2_lowmidlat')
print ('About to read in distribution from ' + datadir + dist_file + ' ...')
dist = np.load(datadir + dist_file)

#if (True): # Debugging purposes
for i in range(nr):
#    i=1 # Debugging purposes
    dist_depth = dist[i]
    cbar_mins[i] = .1
    cbar_maxes[i] = dist_depth.max()
    depth = rr_depth[i]
    # Bin values for this depth
    vrvals, svals = np.meshgrid(vr_bincenters[i], s_bincenters[i])   
    plt.pcolormesh(vrvals, svals, np.transpose(dist_depth), cmap='Greens',\
                   norm=colors.LogNorm(vmin=.1, vmax=dist_depth.max()))
    plt.plot(vr_bincenters[i], np.zeros(nbins_vr), 'k', linewidth=lw)
        # Mark line S' = 0
    plt.xlabel(r'$v_r\ (\rm{m/s})$', fontsize=14)
    plt.ylabel(r'$(S - \langle S\rangle_s)\ (\rm{erg/K/g})$',fontsize=14)
    plt.xlim((np.min(vrvals), 0))
    plt.title(dirname_stripped + r'$\ v_r/S^\prime$' +\
              ' dist (df, lowmidlat), depth = %.2f%%' %(depth*100), fontsize=14)
    # Get ticks everywhere
    plt.minorticks_on()
    plt.tick_params(top='on', right='on', direction='in', which='both')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('counts', rotation=270)
    plt.tight_layout()
    savename = dirname_stripped + '_vr_s_df_depth' +\
        str(i).zfill(3) + '.png'
    print('saving plot at ' + plotdirs[1] + savename + ' ...')
    plt.savefig(plotdirs[1] + savename, dpi=300)
    plt.close()
np.save(datadir + 'vrs_df_cbar_limits_15_to_30.npy', (cbar_maxes, cbar_mins))  
 
# Read in distribution
dist_file = get_widest_range_file(datadir, 'vr_s_dist_df_3_midlat')
print ('About to read in distribution from ' + datadir + dist_file + ' ...')
dist = np.load(datadir + dist_file)

#if (True): # Debugging purposes
for i in range(nr):
#    i=1 # Debugging purposes
    dist_depth = dist[i]
    cbar_mins[i] = .1
    cbar_maxes[i] = dist_depth.max()
    depth = rr_depth[i]
    # Bin values for this depth
    vrvals, svals = np.meshgrid(vr_bincenters[i], s_bincenters[i]) 
    plt.pcolormesh(vrvals, svals, np.transpose(dist_depth), cmap='Greens',\
                   norm=colors.LogNorm(vmin=.1, vmax=dist_depth.max()))
    plt.plot(vr_bincenters[i], np.zeros(nbins_vr), 'k', linewidth=lw)
        # Mark line S' = 0
    plt.xlabel(r'$v_r\ (\rm{m/s})$', fontsize=14)
    plt.ylabel(r'$(S - \langle S\rangle_s)\ (\rm{erg/K/g})$',fontsize=14)
    plt.xlim((np.min(vrvals), 0))
    plt.title(dirname_stripped + r'$\ v_r/S^\prime$' +\
              ' dist (df, midlat), depth = %.2f%%' %(depth*100), fontsize=14)
    # Get ticks everywhere
    plt.minorticks_on()
    plt.tick_params(top='on', right='on', direction='in', which='both')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('counts', rotation=270)
    plt.tight_layout()
    savename = dirname_stripped + '_vr_s_df_depth' +\
        str(i).zfill(3) + '.png'
    print('saving plot at ' + plotdirs[2] + savename + ' ...')
    plt.savefig(plotdirs[2] + savename, dpi=300)
    plt.close()
np.save(datadir + 'vrs_df_cbar_limits_30_to_45.npy', (cbar_maxes, cbar_mins))
   
# Read in distribution
dist_file = get_widest_range_file(datadir, 'vr_s_dist_df_4_midhighlat')
print ('About to read in distribution from ' + datadir + dist_file + ' ...')
dist = np.load(datadir + dist_file)

#if (True): # Debugging purposes
for i in range(nr):
#    i=1 # Debugging purposes
    dist_depth = dist[i]
    cbar_mins[i] = .1
    cbar_maxes[i] = dist_depth.max()
    depth = rr_depth[i]
    # Bin values for this depth
    vrvals, svals = np.meshgrid(vr_bincenters[i], s_bincenters[i])
    plt.pcolormesh(vrvals, svals, np.transpose(dist_depth), cmap='Greens',\
                   norm=colors.LogNorm(vmin=.1, vmax=dist_depth.max()))
    plt.plot(vr_bincenters[i], np.zeros(nbins_vr), 'k', linewidth=lw)
        # Mark line S' = 0
    plt.xlabel(r'$v_r\ (\rm{m/s})$', fontsize=14)
    plt.ylabel(r'$(S - \langle S\rangle_s)\ (\rm{erg/K/g})$',fontsize=14)
    plt.xlim((np.min(vrvals), 0))
    plt.title(dirname_stripped + r'$\ v_r/S^\prime$' +\
              ' dist (df, midhighlat), depth = %.2f%%' %(depth*100), fontsize=14)
    # Get ticks everywhere
    plt.minorticks_on()
    plt.tick_params(top='on', right='on', direction='in', which='both')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('counts', rotation=270)
    plt.tight_layout()
    savename = dirname_stripped + '_vr_s_df_depth' +\
        str(i).zfill(3) + '.png'
    print('saving plot at ' + plotdirs[3] + savename + ' ...')
    plt.savefig(plotdirs[3] + savename, dpi=300)
    plt.close()
np.save(datadir + 'vrs_df_cbar_limits_45_to_60.npy', (cbar_maxes, cbar_mins))

# Read in distribution
dist_file = get_widest_range_file(datadir, 'vr_s_dist_df_5_highlat')
print ('About to read in distribution from ' + datadir + dist_file + ' ...')
dist = np.load(datadir + dist_file)

#if (True): # Debugging purposes
for i in range(nr):
#    i=1 # Debugging purposes
    dist_depth = dist[i]
    cbar_mins[i] = .1
    cbar_maxes[i] = dist_depth.max()
    depth = rr_depth[i]
    # Bin values for this depth
    vrvals, svals = np.meshgrid(vr_bincenters[i], s_bincenters[i])   
    plt.pcolormesh(vrvals, svals, np.transpose(dist_depth), cmap='Greens',\
                   norm=colors.LogNorm(vmin=.1, vmax=dist_depth.max()))
    plt.plot(vr_bincenters[i], np.zeros(nbins_vr), 'k', linewidth=lw)
        # Mark line S' = 0
    plt.xlabel(r'$v_r\ (\rm{m/s})$', fontsize=14)
    plt.ylabel(r'$(S - \langle S\rangle_s)\ (\rm{erg/K/g})$',fontsize=14)
    plt.xlim((np.min(vrvals), 0))
    plt.title(dirname_stripped + r'$\ v_r/S^\prime$' +\
              ' dist (df, highlat), depth = %.2f%%' %(depth*100), fontsize=14)
    # Get ticks everywhere
    plt.minorticks_on()
    plt.tick_params(top='on', right='on', direction='in', which='both')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('counts', rotation=270)
    plt.tight_layout()
    savename = dirname_stripped + '_vr_s_df_depth' +\
        str(i).zfill(3) + '.png'
    print('saving plot at ' + plotdirs[4] + savename + ' ...')
    plt.savefig(plotdirs[4] + savename, dpi=300)
    plt.close()
np.save(datadir + 'vrs_df_cbar_limits_60_to_75.npy', (cbar_maxes, cbar_mins))  
  
# Read in distribution
dist_file = get_widest_range_file(datadir, 'vr_s_dist_df_6_superhighlat')
print ('About to read in distribution from ' + datadir + dist_file + ' ...')
dist = np.load(datadir + dist_file)

#if (True): # Debugging purposes
for i in range(nr):
#    i=1 # Debugging purposes
    dist_depth = dist[i]
    cbar_mins[i] = .1
    cbar_maxes[i] = dist_depth.max()
    depth = rr_depth[i]
    # Bin values for this depth
    vrvals, svals = np.meshgrid(vr_bincenters[i], s_bincenters[i])   
    plt.pcolormesh(vrvals, svals, np.transpose(dist_depth), cmap='Greens',\
                   norm=colors.LogNorm(vmin=.1, vmax=dist_depth.max()))
    plt.plot(vr_bincenters[i], np.zeros(nbins_vr), 'k', linewidth=lw)
        # Mark line S' = 0
    plt.xlabel(r'$v_r\ (\rm{m/s})$', fontsize=14)
    plt.ylabel(r'$(S - \langle S\rangle_s)\ (\rm{erg/K/g})$',fontsize=14)
    plt.xlim((np.min(vrvals), 0))
    plt.title(dirname_stripped + r'$\ v_r/S^\prime$' +\
              ' dist (df, superhighlat), depth = %.2f%%' %(depth*100), fontsize=14)
    # Get ticks everywhere
    plt.minorticks_on()
    plt.tick_params(top='on', right='on', direction='in', which='both')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('counts', rotation=270)
    plt.tight_layout()
    savename = dirname_stripped + '_vr_s_df_depth' +\
        str(i).zfill(3) + '.png'
    print('saving plot at ' + plotdirs[5] + savename + ' ...')
    plt.savefig(plotdirs[5] + savename, dpi=300)
    plt.close()
np.save(datadir + 'vrs_df_cbar_limits_75_to_90.npy', (cbar_maxes, cbar_mins))