import sys, os
import numpy as np
import matplotlib.pyplot as plt
from common import *

# Get the name of the run directory
dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)

# Find the relevant place to store the data, and create the directory if it
# doesn't already exist
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

# Get grid info (if it's not computed already using grid_info.py, this will fail)
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
rr /= 100; ro /= 100; ri /= 100; d /= 100
nr = len(rr)
nt = len(tt)
# Get location of thermal BL
ir_tbl = np.load(datadir + dirname_stripped + '_ir_tbl.npy')
rr_tbl = rr[ir_tbl]
depth_tbl = (ro - rr_tbl)/1.0e6

# Load the trajectories
trajectories = np.load(datadir + dirname_stripped + '_eq_trajectories.npy')
ntraj = len(trajectories)

probs = trajectories[:, -6]
sort_by_probs = np.argsort(probs)
probs = probs[sort_by_probs]
trajectories = trajectories[sort_by_probs, :]

pen_depths = trajectories[:, -1]/1.0e6
pen_depths_log = np.log10(pen_depths)
histvals, binvals, patches = plt.hist(pen_depths_log, weights=probs,\
                                  bins=100, density=True)
pen_depth_log_max_prob = binvals[np.argmax(histvals)]
pen_depth_max_prob = 10**pen_depth_log_max_prob
ax = plt.gca()
xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()
delta_x, delta_y = xmax - xmin, ymax - ymin
yvals = np.linspace(ymin, ymax, 100)
plt.plot(pen_depth_log_max_prob*np.ones(100), yvals, 'k--')
plt.ylim((ymin, ymax))
plt.xlim((xmin, xmax))
ax.text(xmin + 0.1*delta_x, ymax - 0.1*delta_y,\
        'Max pen depth = %.1f Mm' %(np.max(pen_depths)))
ax.text(xmin + 0.1*delta_x, ymax - 0.2*delta_y,\
        'Most prob. = %.1f Mm' %pen_depth_max_prob)
ax.text(xmin + 0.1*delta_x, ymax - 0.3*delta_y,\
        'Depth TBL = %.1f Mm' %depth_tbl)
plt.xlabel('Log10(Penetration depth) (Mm)')
plt.ylabel('Probability density (1/Log10(Mm)')
plt.title(dirname_stripped + ', pen. depth dist.')
plt.tight_layout()
plt.savefig(plotdir + 'traj_hist_pen_depths.png', dpi=300)
plt.close()

# Rossby radii
Rossby_radii = trajectories[:, -4]/1.0e6
histvals, binvals, patches = plt.hist(Rossby_radii, bins=100, weights=probs,\
    density=True)
Rossby_radius_max_prob = binvals[np.argmax(histvals)]
ax = plt.gca()
xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()
delta_x, delta_y = xmax - xmin, ymax - ymin
yvals = np.linspace(ymin, ymax, 100)
plt.plot(Rossby_radius_max_prob*np.ones(100), yvals, 'k--')
ax.text(xmax - 0.5*delta_x, ymax - 0.1*delta_y,\
        'Max Rossby radius = %.1f Mm' %(np.max(Rossby_radii)))
ax.text(xmax - 0.5*delta_x, ymax - 0.2*delta_y,\
        'Most prob. = %.1f Mm' %Rossby_radius_max_prob)
plt.xlim((xmin, xmax))
plt.ylim((ymin, ymax))
#plt.ylim((0, ymax))
plt.xlabel('Rossby radii (Mm)')
plt.ylabel('Probability density (1/Mm)')
plt.title(dirname_stripped + ', Rossby radius dist.')
plt.savefig(plotdir + 'traj_hist_Rossby_radii.png', dpi=300)
plt.close()

# Initial u-values
uvals = trajectories[:, -5]
histvals, binvals, patches = plt.hist(uvals, bins=100, weights=probs,\
                                      density=True)
uval_max_prob = binvals[np.argmax(histvals)]
ax = plt.gca()
xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()
delta_x, delta_y = xmax - xmin, ymax - ymin
yvals = np.linspace(ymin, ymax, 100)
plt.plot(uval_max_prob*np.ones(100), yvals, 'k--')
plt.xlim((xmin, xmax))
plt.ylim((ymin, ymax))
plt.xlabel('u (m/s)')
plt.ylabel('Probability density (s/m)')
plt.title(dirname_stripped + ', u-value dist.')
ax.text(xmax - 0.6*delta_x, ymax - 0.1*delta_y,\
        r'$u = \sqrt{(v_\phi - v_D)^2 + v_r^2}$')
ax.text(xmax - 0.6*delta_x, ymax - 0.2*delta_y,\
        'most prob. = %.1f m/s' %uval_max_prob)
ax.text(xmax - 0.6*delta_x, ymax - 0.3*delta_y,\
        'max(u) = %.1f m/s' %(np.max(uvals)))
plt.savefig(plotdir + 'traj_hist_uvals.png', dpi=300)
plt.close()

# Phase angles
alpha = trajectories[:, -3]/np.pi*180
plt.hist(alpha, weights=probs, bins=100)
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
plt.xlim((-90, 90))
plt.xlabel('Initial Phase angle (arctan((vphi-vD/vr))')
plt.ylabel('Probability density (1/degree)')
plt.title(dirname_stripped + ', phase-angle dist.')
plt.savefig(plotdir + 'traj_hist_phase_angles.png', dpi=300)
plt.close()

# Drift velocities
drift_velocities = trajectories[:, -2]
plt.hist(drift_velocities, weights=probs, bins=100)
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
#plt.xlim((-90, 90))
plt.xlabel(r'$v_D \equiv (S^\prime/c_p)(g/2\Omega_0)\ (m/s)$')
plt.ylabel('Probability density (s/m)')
plt.savefig(plotdir + 'traj_hist_drift_velocities.png', dpi=300)
plt.close()