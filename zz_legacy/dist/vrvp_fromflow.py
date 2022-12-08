import numpy as np
import sys, os
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import colors
import matplotlib.pyplot as plt
from diagnostic_reading import Meridional_Slice
from common import *

# Get the name of the run directory
dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist
datadir = dirname + '/data/'
plotdir = dirname + '/plots/vr_vp_dist/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Get grid info (if it's not computed already using grid_info.py, this will fail)
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr = len(rr)
nt = len(tt)

# Get vavg file
vavg_file = get_widest_range_file(datadir, 'vavg')
vr_av, vt_av, vp_av = np.load(datadir + vavg_file)

# Average vp over low latitudes for consistency
tt_lat = tt*180/np.pi - 90
it1, it2 = np.argmin(np.abs(tt_lat - 15)), np.argmin(np.abs(tt_lat + 15))
vp_av_lowlat = np.mean(vp_av[it1:it2, :], axis=0)

vr2_p,vt2_p,vp2_p,vrvp_p,vrvt_p,vtvp_p,\
    vr2_m,vt2_m, vp2_m, vrvp_m, vrvt_m, vtvp_m, fplus, fminus =\
    np.load(datadir + 'rs_pm.npy')

vr2_t = vr2_p + vr2_m
vp2_t = vp2_p + vp2_m
vrvp_t = vrvp_p + vrvp_m

vr2_lowlat = np.mean(vr2_t[it1:it2, :], axis=0)
vp2_lowlat = np.mean(vp2_t[it1:it2, :], axis=0)
vrvp_lowlat = np.mean(vrvp_t[it1:it2, :], axis=0)

cor_vr_vp_fromflow = vrvp_lowlat/np.sqrt(vr2_lowlat*vp2_lowlat)