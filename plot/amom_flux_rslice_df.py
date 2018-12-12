import numpy as np
import matplotlib.pyplot as plt
import sys, os
from diagnostic_reading import ReferenceState
from common import get_widest_range_file, strip_dirname

dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)

# Read in CLAs (if any) to change default file range to use (by default, this
# will be the last 100 meridional slices)
my_boundstype = 'minmax'
my_min, my_max = (-10, 10) # placeholders
user_specified_minmax = False

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-minmax'):
        my_boundstype = 'manual'
        my_min, my_max = float(args[i+1]), float(args[i+2])
        user_specified_minmax = True
        
# Get data and plot directories, creating the plot directory if it doesn't
# already exist        
datadir = dirname + '/data/'
plotdir = dirname + '/plots/amom_flux/rslice/'

if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)
        
ref = ReferenceState(dirname + '/reference', '')
rho = ref.density

# Get grid info
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr, nt = len(rr), len(tt)

# Get reference state info
rho_2d = rho.reshape((1, nr))
r_2d = rr.reshape((1,nr))
sint_2d = sint.reshape((nt, 1))
rho_r_sint = rho_2d*r_2d*sint_2d

# Read in Reynolds stresses to compute amom_flux
rs_file = get_widest_range_file(datadir, 'rs')
print ('Reading in Reynolds stress from ' + datadir + rs_file + ' ...')
(vr2_p, vt2_p, vp2_p, vrvp_p, vrvt_p, vtvp_p,\
 vr2_m, vt2_m, vp2_m, vrvp_m, vrvt_m, vtvp_m, fplus, fminus) =\
     np.load(datadir + rs_file)

# Also need the average velocities to subtract off the component of the RS
# due to the meridional circulation
vavg_file = get_widest_range_file(datadir, 'vavg')
vr_av, vt_av, vp_av = np.load(datadir + vavg_file)

# Compute the radial and latitudinal angular momentum fluxes in the 
# meridional plane
amom_flux_r_p = rho_r_sint*(vrvp_p - fplus*vr_av*vp_av)
amom_flux_r_m = rho_r_sint*(vrvp_m - fminus*vr_av*vp_av)
amom_flux_r_t = amom_flux_r_p + amom_flux_r_m

# Theta grid in latitude degrees
tt_lat = tt*180/np.pi - 90

# Average over six different latitudes and plot as a function of radius

colors = ['b', 'g', 'r', 'c', 'm', 'y']
xmin, xmax = ri/ro, 1
xvals = np.linspace(xmin, xmax, 100)

# 1, lowlat
it1, it2 = np.argmin(np.abs(tt_lat - 15)), np.argmin(np.abs(tt_lat + 15))
flux_p = np.mean(amom_flux_r_p[it1:it2, :], axis=0)
flux_m = np.mean(amom_flux_r_m[it1:it2, :], axis=0)
flux_t = flux_p + flux_m
flux_df_Busse = flux_p
flux_plumes = flux_m - flux_df_Busse
savename = dirname_stripped + '_amom_flux_r_rslice_plumes_00_to_15.png'
savefile = plotdir + savename
current_color = colors[0]
plt.plot(rr/ro, flux_p, color = current_color, linestyle='--', label='upflow')
plt.plot(rr/ro, flux_m, color = current_color, linestyle=':', label='downflow')
plt.plot(rr/ro, flux_t, color = current_color, label='total')
plt.plot(rr/ro, flux_plumes, color = current_color, linestyle = '-.', label='plumes')
plt.xlim((xmin, xmax))
# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
# Mark the zero line
plt.plot(xvals, np.zeros_like(xvals), 'k')
plt.xlabel(r'$r/r_o$', fontsize=14)
plt.ylabel(r'$\mathcal{F}_r\ (g/s^2)$', fontsize=14)

plt.title(dirname_stripped + ', amom_flux_r, plumes isolated, 00 to 15 degrees lat', fontsize=14)

plt.legend()
plt.tight_layout()
print('Saving figure at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.close()