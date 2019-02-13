# Author: Loren Matilsky
# Created: 02/12/2019
# This script plots the latitudinal forces in the meridional plane,
# treating the advective term differently; sphericity terms from the 
# "v dot grad v" are combined with the Coriolis term in a suggestive 
# way
# To use an AZ_Avgs file
# different than the one assocsiated with the longest averaging range, use
# -usefile [complete name of desired AZ_Avgs file]
# Saves plot in
# [dirname]_theta_force_[first iter]_[last iter].png

import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
sys.path.append(os.environ['pl'])
from rayleigh_diagnostics import ReferenceState
from azavg_util import plot_azav
from common import get_widest_range_file, strip_dirname
from get_parameter import get_parameter
from binormalized_cbar import MidpointNormalize

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)
# Read command-line arguments (CLAs)
my_boundstype = 'manual'
user_specified_minmax = False

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


# See if magnetism is "on"
try:
    magnetism = get_parameter(dirname, 'magnetism')
except:
    magnetism = False # if magnetism wasn't specified, it must be "off"

# Get AZ_Avgs file
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
print ('Getting theta_forces from ' + datadir + AZ_Avgs_file + ' ...')
di = np.load(datadir + AZ_Avgs_file, encoding='latin1').item()

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']

# Get grid info
rr, tt, cost, sint, ri,ro,d = di['rr'], di['tt'], di['cost'],\
        di['sint'], di['ri'], di['ro'], di['d']
nr, nt = len(rr), len(tt)
sint2d = sint.reshape((nt, 1))
cost2d = cost.reshape((nt, 1))
rr2d = rr.reshape((1, nr))

# Get density
ref = ReferenceState(dirname + '/reference')
rho = ref.density
rho2d = rho.reshape((1, nr))

# Compute "built-in" forces first
ind_prs = lut[1238]
ind_visc = lut[1229]
theta_force_prs = vals[:, :, ind_prs]
theta_force_visc = vals[:, :, ind_visc]

# Now "custom" inertial forces
ind_vgrad_vt = lut[2212]
ind_pvgrad_pvt = lut[2213]
ind_vrvt = lut[2206]
ind_pvrpvt = lut[2209]
ind_vp = lut[3]
ind_phivsq = lut[424]

vgrad_vt = vals[:, :, ind_vgrad_vt]
pvgrad_pvt = vals[:, :, ind_pvgrad_pvt]
mvgrad_mvt = vgrad_vt - pvgrad_pvt
vrvt = vals[:, :, ind_vrvt]
pvrpvt = vals[:, :, ind_pvrpvt]
mvrmvt = vrvt - pvrpvt
vp = vals[:, :, ind_vp]
phivsq = vals[:, :, ind_phivsq]

vgradv_term_mean = -rho2d*mvgrad_mvt
vgradv_term_prime = -rho2d*pvgrad_pvt
vrvt_term_mean = -rho2d*mvrmvt/rr2d
vrvt_term_prime = -rho2d*pvrpvt/rr2d
phivsq_term = cost2d/sint2d/rr2d*rho2d*phivsq
om0 = get_parameter(dirname, 'angular_velocity')
om = om0 + vp/rr2d/sint2d
diffrot_term = rho2d*sint2d*cost2d*rr2d*(om**2 - om0**2)

vgradv_tot = vgradv_term_mean + vgradv_term_prime + vrvt_term_mean +\
        vrvt_term_prime + phivsq_term + diffrot_term

theta_force_tot = vgradv_tot + theta_force_prs + theta_force_visc

max_sig = max(np.std(vgradv_term_mean), np.std(vgradv_term_prime),\
            np.std(vrvt_term_mean), np.std(vrvt_term_prime), 
            np.std(phivsq_term), np.std(diffrot_term), np.std(vgradv_tot),
              np.std(theta_force_prs), np.std(theta_force_visc))

if magnetism:
    ind_mag = lut[1249]
    theta_force_mag = vals[:, :, ind_mag]       
    max_sig = max(max_sig, np.std(theta_force_mag))
    theta_force_tot += theta_force_mag
    
if not user_specified_minmax: 
    my_min, my_max = -3*max_sig, 3*max_sig

# Set up the actual figure from scratch
fig_width_inches = 7 # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1/8 # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_top_inches = 2 # wider top margin to accommodate subplot titles AND metadata
margin_subplot_top_inches = 1 # margin to accommodate just subplot titles
nplots = 10 + magnetism
ncol = 3 # put three plots per row
nrow = np.int(np.ceil(nplots/3))

subplot_width_inches = (fig_width_inches - (ncol + 1)*margin_inches)/ncol
    # Make the subplot width so that ncol subplots fit together side-by-side
    # with margins in between them and at the left and right.
subplot_height_inches = 2*subplot_width_inches # Each subplot should have an
    # aspect ratio of y/x = 2/1 to accommodate meridional planes. 
fig_height_inches = nrow*subplot_height_inches + margin_top_inches +\
    (nrow - 1)*margin_subplot_top_inches + margin_inches 
    # Room for titles on each row and a regular margin on the bottom
fig_aspect = fig_height_inches/fig_width_inches

# "Margin" in "figure units"; figure units extend from 0 to 1 in BOTH 
# directions, so unitless dimensions of margin will be different in x and y
# to force an equal physical margin
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_subplot_top = margin_subplot_top_inches/fig_height_inches

# Subplot dimensions in figure units
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches


theta_forces = [vgradv_term_mean, vgradv_term_prime, vrvt_term_mean,\
        vrvt_term_prime, phivsq_term, diffrot_term, vgradv_tot,\
        theta_force_prs, theta_force_visc, theta_force_tot]

titles =\
[r'$-\overline{\rho}\langle \mathbf{v}\rangle\cdot\nabla\langle v_\theta\langle$',\
r'$-\overline{\rho}\langle\mathbf{v}^\prime \cdot\nabla v_\theta^\prime\rangle$',\
r'$-\overline{\rho}\langle v_r\rangle\langle v_\theta\rangle/r$',\
r'$-\overline{\rho}\langle v_r^\prime v_\theta^\prime\rangle/r$',\
r'$(\cot{\theta}/r)\overline{\rho}\langle(v_\phi^\prime)^2\rangle$',\
r'$(1/2)\overline{\rho}r\sin{2\theta}(\Omega^2-\Omega_0^2)$',\
 r'$(\mathbf{f}_{\rm{adv,\ tot}})_\theta$',\
 r'$(\mathbf{f}_{\rm{p}})_\theta$',\
 r'$(\mathbf{f}_{\rm{v}})_\theta$',\
 r'$(\mathbf{f}_{\rm{tot}})_\theta$']

units = r'$\rm{g}\ \rm{cm}^{-2}\ \rm{s}^{-2}$'

if magnetism:
    theta_forces.insert(4, theta_force_mag)
    titles.insert(4, r'$(\mathbf{f}_{\rm{mag}})_\theta$')

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

for iplot in range(nplots):
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1 - margin_top - subplot_height - \
            (iplot//ncol)*(subplot_height + margin_subplot_top)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    plot_azav (fig, ax, theta_forces[iplot], rr, cost, sint,\
           units = units,\
           boundstype = my_boundstype, caller_minmax = (my_min, my_max),\
           norm=MidpointNormalize(0))

    ax.set_title(titles[iplot], verticalalignment='bottom', **csfont)

# Put some metadata in upper left
fsize = 12
fig.text(margin_x, 1 - 0.1*margin_top, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.3*margin_top, 'Latitudinal force balance (zonally averaged)',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.5*margin_top,\
         str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8),\
         ha='left', va='top', fontsize=fsize, **csfont)

savefile = plotdir + dirname_stripped + '_theta_forces-2_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

print ('Saving theta_forces (more terms) at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.show()
