import numpy as np
import matplotlib.pyplot as plt
from common import get_widest_range_file, strip_dirname
from diagnostic_reading import ReferenceState
import os, sys

dirname = sys.argv[1]
datadir = dirname + '/data/'
dirname_stripped = strip_dirname(dirname)

# density
ref = ReferenceState('reference', dirname + '/')
rho = ref.density

args = sys.argv[2:]
nargs = len(args)

# Grid info:
rr, tt, cost, sint, rr_depth, ri, ro, d = np.load(datadir + 'grid_info.npy')
nr = len(rr)
nt = len(tt)
rrn = rr/ro

lat_range = 1
rrn0 = 0.976 # corresponds to 10% depth

my_min, my_max = None, None
for i in range(nargs):
    arg = args[i]
    if arg == '-lat':
        lat_range = int(args[i+1])
    elif arg == '-r':
        rrn0 = float(args[i+1])
    elif arg == '-minmax':
        my_min, my_max = float(args[i+1]), float(args[i+2])
        
ir0 = np.argmin(np.abs(rrn - rrn0))

latrange_dict = {1:'vr_vp_dist_full_1_lowlat', 2:'vr_vp_dist_full_2_lowmidlat', 3:'vr_vp_dist_full_3_midlat',\
                 4:'vr_vp_dist_full_4_midhighlat', 5:'vr_vp_dist_full_5_highlat', 6:'vr_vp_dist_full_6_superhighlat'}
latvals_dict = {1:(0, 15), 2:(15, 30), 3:(30, 45), 4:(45, 60), 5:(60, 75), 6:(75, 90)}

dist_file = get_widest_range_file(datadir, latrange_dict[lat_range])
print ('Reading %s ...' %dist_file)
dist = np.load(datadir + dist_file)

# Get the bin structure for the distribution
vr_bincenters, vp_bincenters = np.load(datadir + 'vrvp_bincenters.npy')
nbins_vr, nbins_vp = len(vr_bincenters[0]), len(vp_bincenters[0])

# Get vavg file and average over the appropriate latitude region
vavg_file = get_widest_range_file(datadir, 'vavg')
vr_av, vt_av, vp_av = np.load(datadir + vavg_file)
vr_av /= 100; vt_av /= 100; vp_av /= 100 # Get average velocities in m/s

lat_low, lat_high = latvals_dict[lat_range]
tt_lat = 90 - tt*180/np.pi
it1, it2 = np.argmin(np.abs(tt_lat + lat_high)), np.argmin(np.abs(tt_lat + lat_low))
it3, it4 = np.argmin(np.abs(tt_lat - lat_low)), np.argmin(np.abs(tt_lat - lat_high))

print ('looking at lat range it1, it2, it3, it4 = ', it1, it2, it3, it4)
print (' and depth rr/ro = ', rrn[ir0], ' ir0 = ', ir0)

# DONT average first...doesn't make sense!
#vr_av0 = 0.5*(np.mean(vr_av[it1:it2+1, ir0]) + np.mean(vr_av[it3:it4+1, ir0]))
#vp_av0 = 0.5*(np.mean(vp_av[it1:it2+1, ir0]) + np.mean(vp_av[it3:it4+1, ir0]))
vr_av0 = vr_av
print ('vr_av0, vp_av0 = ', vr_av0, vp_av0)

# Get local distribution
dist0 = dist[ir0]
vrvals0 = vr_bincenters[ir0]
vpvals0 = vp_bincenters[ir0]

vrvals0_2d, vpvals0_2d = np.meshgrid(vrvals0, vpvals0, indexing='ij')

# Get velocity fluctuations
vr_fluc0 = vrvals0_2d - vr_av0
vp_fluc0 = vpvals0_2d - vp_av0

# Compute cdf's of amom flux
amom0_dist = vp_fluc0*vr_fluc0*dist0/np.sum(dist0)*rho[ir0]*rr[ir0]*1e4
amom0_dist_r = np.sum(amom0_dist, axis=1)
nvr0 = len(vrvals0)

where_up = np.where(vrvals0 > 0)
where_down = np.where(vrvals0 < 0)

amom0_dist_r_up = amom0_dist_r[where_up]
amom0_dist_r_down = amom0_dist_r[where_down]

vrvals0_up = vrvals0[where_up]
vrvals0_down = vrvals0[where_down]
nvr_up, nvr_down = len(vrvals0_up), len(vrvals0_down)

amom0_cdf_up = np.zeros(nvr_up)
for i in range(nvr_up):
    amom0_cdf_up[i] = np.sum(amom0_dist_r_up[:i])

amom0_cdf_down = np.zeros(nvr_down)
for i in range(nvr_down):
    amom0_cdf_down[i] = np.sum(amom0_dist_r_down[i:])

plt.plot(vrvals0_up, amom0_cdf_up, 'r')
plt.plot(vrvals0_down, amom0_cdf_down, 'b')

datamin, datamax = min(np.min(amom0_cdf_up), np.min(amom0_cdf_down)), max(np.max(amom0_cdf_up), np.max(amom0_cdf_down))
diff = datamax - datamin
buffer = 0.1*diff
ymin, ymax = datamin - buffer, datamax + buffer
ymin = min(ymin, 0)
if not my_min is None:
    ymin, ymax = my_min, my_max
yvals = np.linspace(ymin, ymax, 100)
xmin, xmax = np.min(vrvals0), np.max(vrvals0)
xvals = np.linspace(xmin, xmax, 100)
plt.xlim((xmin, xmax))
plt.ylim((ymin, ymax))
plt.plot(xvals, np.zeros(100), 'k', linewidth=0.5)
plt.plot(np.zeros(100), yvals, 'k', linewidth=0.5)
plt.xlabel('vr (m/s)')
plt.ylabel("cdf(F_r) (g/s^2)")
plt.title('ir = %i, r/ro = %0.3f, (lat_low, lat_high) = (%i, %i)' %(ir0, rrn[ir0], lat_low, lat_high))
plt.minorticks_on()
plt.tick_params(top=True, left=True, right=True, bottom=True, direction='in', which='both')
plt.tight_layout()

plotname = dirname_stripped + '_amomflux_cdf_' + ('latrange%i_' %lat_range) + ('rrn%0.3f' %rrn[ir0]) + '.png'

savedir = '/altair/loma3853/rayleigh/all_plots/amomflux_cdf/' + dirname_stripped + '/'
if not os.path.isdir(savedir):
    os.makedirs(savedir)
plt.savefig(savedir + plotname, dpi=300)
