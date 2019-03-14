import numpy as np
import pickle
import matplotlib.pyplot as plt
import sys
import os
from azavg_util import plot_azav
from binormalized_cbar import MidpointNormalize
from diagnostic_reading import ReferenceState
from common import get_widest_range_file, strip_dirname

dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)

datadir = dirname + '/data/'
plotdir = dirname + '/plots/amom_flux/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Read in CLAs (if any) to change default file range to use (by default, this
# will be the last 100 meridional slices)
my_boundstype = 'minmax'
my_min, my_max = (-10, 10) # placeholders
user_specified_minmax = False
showplot = False

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-minmax'):
        my_boundstype = 'manual'
        my_min, my_max = float(args[i+1]), float(args[i+2])
        user_specified_minmax = True
    if (arg == '-show'):
        showplot = True
    if (arg == '-nlevs'):
        my_nlevs = int(args[i+1])
        
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
print ('Getting Reynolds stress from ' + datadir + rs_file + ' ...')
(vr2_p, vt2_p, vp2_p, vrvp_p, vrvt_p, vtvp_p,\
     vr2_m, vt2_m, vp2_m, vrvp_m, vrvt_m, vtvp_m, fplus, fminus) =\
     np.load(datadir + rs_file)

# Also need the average velocities to subtract off the component of the RS
# due to the meridional circulation
vavg_file = get_widest_range_file(datadir, 'vavg')
vr_av, vt_av, vp_av = np.load(datadir + vavg_file)
    
amom_flux_t_p = rho_r_sint*(vtvp_p - fplus*vt_av*vp_av)
amom_flux_t_m = rho_r_sint*(vtvp_m - fminus*vt_av*vp_av)
amom_flux_t_t = amom_flux_t_p + amom_flux_t_m

# Plot latitudinal amom_flux, upflows
fig, ax = plt.subplots()
plot_azav (fig, ax, amom_flux_t_p, rr, cost, sint, units = r'$g/s^2$',\
           plotcontours=False, boundstype = my_boundstype,\
           caller_minmax = (my_min, my_max), norm=MidpointNormalize(0))
plt.title(dirname_stripped + r'$,\ \mathcal{F}_{\theta,\ uf}$')
plt.tight_layout()
# Save the plot
savefile = plotdir + dirname_stripped + '_amom_flux_t_uf.png'
print ('Saving figure at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.close()

# Plot latitudinal amom_flux, downflows
fig, ax = plt.subplots()
plot_azav (fig, ax, amom_flux_t_m, rr, cost, sint, units = r'$g/s^2$',\
           plotcontours=False, boundstype = my_boundstype,\
           caller_minmax = (my_min, my_max), norm=MidpointNormalize(0))
plt.title(dirname_stripped + r'$,\ \mathcal{F}_{\theta,\ df}$')
plt.tight_layout()
# Save the plot
savefile = plotdir + dirname_stripped + '_amom_flux_t_df.png'
print ('Saving figure at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.close()

# Plot latitudinal amom_flux, total
fig, ax = plt.subplots()
plot_azav (fig, ax, amom_flux_t_t, rr, cost, sint, units = r'$g/s^2$',\
           plotcontours=False, boundstype = my_boundstype,\
           caller_minmax = (my_min, my_max), norm=MidpointNormalize(0))
plt.title(dirname_stripped + r'$,\ \mathcal{F}_{\theta,\ tot}$')
plt.tight_layout()
# Save the plot
savefile = plotdir + dirname_stripped + '_amom_flux_t_tot.png'
print ('Saving figure at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.close()
