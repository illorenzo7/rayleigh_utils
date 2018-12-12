import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from azavg_util import plot_azav
from binormalized_cbar import MidpointNormalize
from diagnostic_reading import ReferenceState
from common import get_widest_range_file, strip_dirname

# Read in directory of the Rayleigh run
dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)

# Read in CLAs (if any) to change default file range to use (by default, this
# will be the last 100 meridional slices)
my_boundstype = 'minmax'
my_min, my_max = (-10, 10) # placeholders
user_specified_minmax = False
fromdist = False

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-minmax'):
        my_boundstype = 'manual'
        my_min, my_max = float(args[i+1]), float(args[i+2])
        user_specified_minmax = True
    if (arg == '-fromdist'):
        fromdist = True

datadir = dirname + '/data/'
plotdir = dirname + '/plots/amom_flux/'
if (fromdist):
    plotdir += 'fromdist/'   
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
if (not fromdist):
    rs_file = get_widest_range_file(datadir, 'rs')
    print ('Reading in Reynolds stress from ' + datadir + rs_file + ' ...')
    (vr2_p, vt2_p, vp2_p, vrvp_p, vrvt_p, vtvp_p,\
     vr2_m, vt2_m, vp2_m, vrvp_m, vrvt_m, vtvp_m, fplus, fminus) =\
         np.load(datadir + rs_file)
else:
    rs_file = get_widest_range_file(datadir, 'rs_fromdist')
    print ('Reading in Reynolds stress from ' + datadir + rs_file + ' ...')
    vr2_p, vp2_p, vrvp_p, vr2_m, vp2_m, vrvp_m, fplus, fminus =\
        np.load(datadir + rs_file)


# Also need the average velocities to subtract off the component of the RS
# due to the meridional circulation
vavg_file = get_widest_range_file(datadir, 'vavg')
vr_av, vt_av, vp_av = np.load(datadir + vavg_file)
    
amom_flux_r_p = rho_r_sint*(vrvp_p - fplus*vr_av*vp_av)
amom_flux_r_m = rho_r_sint*(vrvp_m - fminus*vr_av*vp_av)
amom_flux_r_t = amom_flux_r_p + amom_flux_r_m

# Plot radial amom_flux, upflows
fig, ax = plt.subplots()
plot_azav (fig, ax, amom_flux_r_p, rr, cost, sint, units = r'$g/s^2$',\
           plotcontours=False, boundstype = my_boundstype,\
           caller_minmax = (my_min, my_max), norm=MidpointNormalize(0))
if (fromdist):
    plt.title(dirname_stripped + r'$,\ \mathcal{F}_{r,\ uf}$' + ' from dist')
else:
    plt.title(dirname_stripped + r'$,\ \mathcal{F}_{r,\ uf}$')
plt.tight_layout()
# Save the plot
if (fromdist):
    savefile = plotdir + dirname_stripped + '_amom_flux_r_uf_fromdist.png'
else:
    savefile = plotdir + dirname_stripped + '_amom_flux_r_uf.png'
print ('Saving figure at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.close()

# Plot radial amom_flux, downflows
fig, ax = plt.subplots()
plot_azav (fig, ax, amom_flux_r_m, rr, cost, sint, units = r'$g/s^2$',\
           plotcontours=False, boundstype = my_boundstype,\
           caller_minmax = (my_min, my_max), norm=MidpointNormalize(0))
if (fromdist):
    plt.title(dirname_stripped + r'$,\ \mathcal{F}_{r,\ df}$' + ' from dist')
else:
    plt.title(dirname_stripped + r'$,\ \mathcal{F}_{r,\ df}$')
plt.tight_layout()
# Save the plot
if (fromdist):
    savefile = plotdir + dirname_stripped + '_amom_flux_r_df_fromdist.png'
else:
    savefile = plotdir + dirname_stripped + '_amom_flux_r_df.png'
print ('Saving figure at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.close()

# Plot radial amom_flux, total
fig, ax = plt.subplots()
plot_azav (fig, ax, amom_flux_r_t, rr, cost, sint, units = r'$g/s^2$',\
           plotcontours=False, boundstype = my_boundstype,\
           caller_minmax = (my_min, my_max), norm=MidpointNormalize(0))
if (fromdist):
    plt.title(dirname_stripped + r'$,\ \mathcal{F}_{r,\ tot}$' + ' from dist')
else:
    plt.title(dirname_stripped + r'$,\ \mathcal{F}_{r,\ tot}$')
plt.tight_layout()
# Save the plot
if (fromdist):
    savefile = plotdir + dirname_stripped + '_amom_flux_r_tot_fromdist.png'
else:
    savefile = plotdir + dirname_stripped + '_amom_flux_r_tot.png'
print ('Saving figure at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.close()