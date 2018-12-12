import numpy as np
from azavg_util import plot_azav
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys, os
from binormalized_cbar import MidpointNormalize
from common import get_widest_range_file, get_iters_from_file, strip_dirname

dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

user_specified_minmax = False
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        user_min = float(args[i+1])
        user_max = float(args[i+2])
        user_specified_minmax = True

datadir = dirname + '/data/'
plotdir = dirname + '/plots/'

if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Get grid info
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr, nt = len(rr), len(tt)

eflux_merplane_file = get_widest_range_file(datadir, 'eflux_radial_merplane')
vflux, kflux, cflux, hflux, eflux, tflux = np.load(datadir + eflux_merplane_file)

dphi, dt, dr = np.load(datadir + 'differentials.npy')
rsint = rr.reshape((1, nr))*sint.reshape((nt, 1))
eflux_av = np.sum(eflux*dr.reshape((1,nr)), axis=1)/np.sum(dr)
it_half = nt//2 
eflux_av_90 = 0.5*(eflux_av[:it_half] + eflux_av[it_half:][::-1])

savename = plotdir + dirname_stripped + '_enth_flux_vs_latitude.png'

tt_lat = tt*180/np.pi - 90
it1, it2 = np.argmin(np.abs(tt_lat - 75)), np.argmin(np.abs(tt_lat + 75))
latmax = float(sys.argv[2])
it_max = np.argmin(np.abs(tt_lat - latmax))


plt.plot(tt_lat[it_max:it_half], eflux_av_90[it_max:])
plt.ylim((0, 1.2*np.max(eflux_av_90[it_max:])))

plt.savefig(savename, dpi=300)
plt.close()
