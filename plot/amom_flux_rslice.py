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
        
# Get data and plot directories, creating the plot directory if it doesn't
# already exist        
datadir = dirname + '/data/'
plotdir = dirname + '/plots/amom_flux/rslice/'
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


# Compute the radial and latitudinal angular momentum fluxes in the 
# meridional plane
amom_flux_r_p = rho_r_sint*vrvp_p # DO NOT subtract off meridional circulation parts!
 											# This was already done in rs.py/rs_fromdist.py -LIM, 07/22/2018
amom_flux_r_m = rho_r_sint*vrvp_m
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
savename = dirname_stripped + '_amom_flux_r_rslice_00_to_15.png'
if (fromdist):
    savename = dirname_stripped + '_amom_flux_r_rslice_00_to_15_fromdist.png'
savefile = plotdir + savename
current_color = colors[0]
plt.plot(rr/ro, flux_p, color = current_color, linestyle='--', label='upflow')
plt.plot(rr/ro, flux_m, color = current_color, linestyle=':', label='downflow')
plt.plot(rr/ro, flux_t, color = current_color, label='total')
plt.xlim((xmin, xmax))
plt.ylim((-1.5e16, 2e16))
# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
# Mark the zero line
plt.plot(xvals, np.zeros_like(xvals), 'k')
plt.xlabel(r'$r/r_o$', fontsize=14)
plt.ylabel(r'$\mathcal{F}_r\ (g/s^2)$', fontsize=14)

if (not fromdist):
    plt.title(dirname_stripped + ', amom_flux_r, 00 to 15 degrees lat', fontsize=14)
else:
    plt.title(dirname_stripped + ', amom_flux_r, from dist, 00 to 15 degrees lat',\
              fontsize=14)
plt.legend()
plt.tight_layout()
print('Saving figure at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.close()

# 2, lowmidlat
it1, it2 = np.argmin(np.abs(tt_lat - 30)), np.argmin(np.abs(tt_lat - 15))
it3, it4 = np.argmin(np.abs(tt_lat + 15)), np.argmin(np.abs(tt_lat + 30))
flux_p = 0.5*(np.mean(amom_flux_r_p[it1:it2, :], axis=0) +\
              np.mean(amom_flux_r_p[it3:it4, :], axis=0))
flux_m = 0.5*(np.mean(amom_flux_r_m[it1:it2, :], axis=0) +\
              np.mean(amom_flux_r_m[it3:it4, :], axis=0))
flux_t = flux_p + flux_m
savename = dirname_stripped + '_amom_flux_r_rslice_15_to_30.png'
if (fromdist):
    savename = dirname_stripped + '_amom_flux_r_rslice_15_to_30_fromdist.png'
savefile = plotdir + savename
current_color = colors[1]
plt.plot(rr/ro, flux_p, color = current_color, linestyle='--', label='upflow')
plt.plot(rr/ro, flux_m, color = current_color, linestyle=':', label='downflow')
plt.plot(rr/ro, flux_t, color = current_color, label='total')
plt.xlim((xmin, xmax))
# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
# Mark the zero line
plt.plot(xvals, np.zeros_like(xvals), 'k')
plt.xlabel(r'$r/r_o$', fontsize=14)
plt.ylabel(r'$\mathcal{F}_r\ (g/s^2)$', fontsize=14)
if (not fromdist):
    plt.title(dirname_stripped + ', amom_flux_r, 15 to 30 degrees lat', fontsize=14)
else:
    plt.title(dirname_stripped + ', amom_flux_r, from dist, 15 to 30 degrees lat',\
              fontsize=14)
plt.legend()
plt.tight_layout()
print('Saving figure at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.close()


# 3, midlat
it1, it2 = np.argmin(np.abs(tt_lat - 45)), np.argmin(np.abs(tt_lat - 30))
it3, it4 = np.argmin(np.abs(tt_lat + 30)), np.argmin(np.abs(tt_lat + 45))
flux_p = 0.5*(np.mean(amom_flux_r_p[it1:it2, :], axis=0) +\
              np.mean(amom_flux_r_p[it3:it4, :], axis=0))
flux_m = 0.5*(np.mean(amom_flux_r_m[it1:it2, :], axis=0) +\
              np.mean(amom_flux_r_m[it3:it4, :], axis=0))
flux_t = flux_p + flux_m
savename = dirname_stripped + '_amom_flux_r_rslice_30_to_45.png'
if (fromdist):
    savename = dirname_stripped + '_amom_flux_r_rslice_30_to_45_fromdist.png'
savefile = plotdir + savename
current_color = colors[2]
plt.plot(rr/ro, flux_p, color = current_color, linestyle='--', label='upflow')
plt.plot(rr/ro, flux_m, color = current_color, linestyle=':', label='downflow')
plt.plot(rr/ro, flux_t, color = current_color, label='total')
plt.xlim((xmin, xmax))
# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
# Mark the zero line
plt.plot(xvals, np.zeros_like(xvals), 'k')
plt.xlabel(r'$r/r_o$', fontsize=14)
plt.ylabel(r'$\mathcal{F}_r\ (g/s^2)$', fontsize=14)
if (not fromdist):
    plt.title(dirname_stripped + ', amom_flux_r, 30 to 45 degrees lat', fontsize=14)
else:
    plt.title(dirname_stripped + ', amom_flux_r, from dist, 30 to 45 degrees lat',\
              fontsize=14)
plt.legend()
plt.tight_layout()
print('Saving figure at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.close()


# 4, midhighlat
it1, it2 = np.argmin(np.abs(tt_lat - 60)), np.argmin(np.abs(tt_lat - 45))
it3, it4 = np.argmin(np.abs(tt_lat + 45)), np.argmin(np.abs(tt_lat + 60))
flux_p = 0.5*(np.mean(amom_flux_r_p[it1:it2, :], axis=0) +\
              np.mean(amom_flux_r_p[it3:it4, :], axis=0))
flux_m = 0.5*(np.mean(amom_flux_r_m[it1:it2, :], axis=0) +\
              np.mean(amom_flux_r_m[it3:it4, :], axis=0))
flux_t = flux_p + flux_m
savename = dirname_stripped + '_amom_flux_r_rslice_45_to_60.png'
if (fromdist):
    savename = dirname_stripped + '_amom_flux_r_rslice_45_to_60_fromdist.png'
savefile = plotdir + savename
current_color = colors[3]
plt.plot(rr/ro, flux_p, color = current_color, linestyle='--', label='upflow')
plt.plot(rr/ro, flux_m, color = current_color, linestyle=':', label='downflow')
plt.plot(rr/ro, flux_t, color = current_color, label='total')
plt.xlim((xmin, xmax))
plt.ylim((-1.2e16, 0.2e16))
# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
# Mark the zero line
plt.plot(xvals, np.zeros_like(xvals), 'k')
plt.xlabel(r'$r/r_o$', fontsize=14)
plt.ylabel(r'$\mathcal{F}_r\ (g/s^2)$', fontsize=14)
if (not fromdist):
    plt.title(dirname_stripped + ', amom_flux_r, 45 to 60 degrees lat', fontsize=14)
else:
    plt.title(dirname_stripped + ', amom_flux_r, from dist, 45 to 60 degrees lat',\
              fontsize=14)
plt.legend()
plt.tight_layout()
print('Saving figure at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.close()

# 5, highlat
it1, it2 = np.argmin(np.abs(tt_lat - 75)), np.argmin(np.abs(tt_lat - 60))
it3, it4 = np.argmin(np.abs(tt_lat + 60)), np.argmin(np.abs(tt_lat + 75))
flux_p = 0.5*(np.mean(amom_flux_r_p[it1:it2, :], axis=0) +\
              np.mean(amom_flux_r_p[it3:it4, :], axis=0))
flux_m = 0.5*(np.mean(amom_flux_r_m[it1:it2, :], axis=0) +\
              np.mean(amom_flux_r_m[it3:it4, :], axis=0))
flux_t = flux_p + flux_m
savename = dirname_stripped + '_amom_flux_r_rslice_60_to_75.png'
if (fromdist):
    savename = dirname_stripped + '_amom_flux_r_rslice_60_to_75_fromdist.png'
savefile = plotdir + savename
current_color = colors[4]
plt.plot(rr/ro, flux_p, color = current_color, linestyle='--', label='upflow')
plt.plot(rr/ro, flux_m, color = current_color, linestyle=':', label='downflow')
plt.plot(rr/ro, flux_t, color = current_color, label='total')
plt.xlim((xmin, xmax))
# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
# Mark the zero line
plt.plot(xvals, np.zeros_like(xvals), 'k')
plt.xlabel(r'$r/r_o$', fontsize=14)
plt.ylabel(r'$\mathcal{F}_r\ (g/s^2)$', fontsize=14)
if (not fromdist):
    plt.title(dirname_stripped + ', amom_flux_r, 60 to 75 degrees lat', fontsize=14)
else:
    plt.title(dirname_stripped + ', amom_flux_r, from dist, 60 to 75 degrees lat',\
              fontsize=14)
plt.legend()
plt.tight_layout()
print('Saving figure at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.close()

# 6, superhighlat
it1, it2 = np.argmin(np.abs(tt_lat - 90)), np.argmin(np.abs(tt_lat - 75))
it3, it4 = np.argmin(np.abs(tt_lat + 75)), np.argmin(np.abs(tt_lat + 90))
flux_p = 0.5*(np.mean(amom_flux_r_p[it1:it2, :], axis=0) +\
              np.mean(amom_flux_r_p[it3:it4, :], axis=0))
flux_m = 0.5*(np.mean(amom_flux_r_m[it1:it2, :], axis=0) +\
              np.mean(amom_flux_r_m[it3:it4, :], axis=0))
flux_t = flux_p + flux_m
savename = dirname_stripped + '_amom_flux_r_rslice_75_to_90.png'
if (fromdist):
    savename = dirname_stripped + '_amom_flux_r_rslice_75_to_90_fromdist.png'
savefile = plotdir + savename
current_color = colors[5]
plt.plot(rr/ro, flux_p, color = current_color, linestyle='--', label='upflow')
plt.plot(rr/ro, flux_m, color = current_color, linestyle=':', label='downflow')
plt.plot(rr/ro, flux_t, color = current_color, label='total')
plt.xlim((xmin, xmax))
# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
# Mark the zero line
plt.plot(xvals, np.zeros_like(xvals), 'k')
plt.xlabel(r'$r/r_o$', fontsize=14)
plt.ylabel(r'$\mathcal{F}_r\ (g/s^2)$', fontsize=14)
if (not fromdist):
    plt.title(dirname_stripped + ', amom_flux_r, 75 to 90 degrees lat', fontsize=14)
else:
    plt.title(dirname_stripped + ', amom_flux_r, from dist, 75 to 90 degrees lat',\
              fontsize=14)
plt.legend()
plt.tight_layout()
print('Saving figure at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.close()
