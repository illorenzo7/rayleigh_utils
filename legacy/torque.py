# LEGACY TORQUE COMPUTATION (only good for n3 cases, n4, n5, etc.)
# NEEDS UPDATED DOCS AND TESTING
from derivs import drad, dth
import numpy as np
import sys, os
from diagnostic_reading import ReferenceState
from common import get_widest_range_file, get_iters_from_file, strip_dirname
from get_parameter import get_parameter

dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

datadir = dirname + '/data/'
plotdir = dirname + '/plots/'

if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Read command-line arguments (CLAs)
my_boundstype = 'minmax'
my_min, my_max = -10, 10
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
        

ref = ReferenceState(dirname + '/reference', '')
rho = ref.density

# Get grid info
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr, nt = len(rr), len(tt)

rho_2d = rho.reshape((1, nr))
r_2d = rr.reshape((1,nr))
sint_2d = sint.reshape((nt, 1))
rho_r_sint = rho_2d*r_2d*sint_2d

# Get the average velocity to calculate 
# meridional circulation torque and viscous torque
vavg_file = get_widest_range_file(datadir, 'vavg')
iter1, iter2 = get_iters_from_file(vavg_file)
vr_av, vt_av, vp_av = np.load(datadir + vavg_file)

## Get Reynolds stresses for Reynolds stress torque
#rs_file = get_widest_range_file(datadir, 'rs')
#vr2_p,vt2_p,vp2_p,vrvp_p,vrvt_p,vtvp_p,\
#    vr2_m,vt2_m,vp2_m, vrvp_m, vrvt_m, vtvp_m, fplus, fminus =\
#        np.load(datadir + rs_file)
#vrvp = vrvp_m + vrvp_p
#vtvp = vtvp_m + vtvp_p
#vrvp -= vr_av*vp_av
#vtvp -= vt_av*vp_av

# Torque due to meridional circulation
# First, fluxes
om0 = get_parameter(dirname, 'angular_velocity')
L = r_2d*sint_2d*(om0*r_2d*sint_2d + vp_av)
torque_mc_r = -rho_2d*vr_av*drad(L, rr)
torque_mc_t = -rho_2d*vt_av*dth(L, tt)/r_2d
torque_mc = torque_mc_r + torque_mc_t

# alternatively:
#f_mc_r = rho_r_sint*vr_av*(vp_av + om0*r_2d*sint_2d)
#f_mc_t = rho_r_sint*vt_av*(vp_av + om0*r_2d*sint_2d)
#torque_mc3 = -(1./r_2d**2*drad(r_2d**2*f_mc_r,rr) +\
#        1./r_2d/sint_2d*dth(sint_2d*f_mc_t,tt))

# Torque due to Reynolds stress
f_rs_file = get_widest_range_file(datadir, 'amomflux_rs')
f_rs_r, f_rs_t = np.load(datadir + f_rs_file)
torque_rs = -(1./r_2d**2*drad(r_2d**2*f_rs_r,rr) +\
        1./r_2d/sint_2d*dth(sint_2d*f_rs_t,tt))

# Torque due to viscosity
nu = get_parameter(dirname, 'nu_top')
f_visc_t = dth(vp_av/sint_2d, tt)
f_visc_t = -nu*sint_2d**2*rho_2d*f_visc_t

f_visc_r = drad(vp_av/r_2d, rr)
f_visc_r = -nu*r_2d**2*sint_2d*rho_2d*f_visc_r

torque_visc = -(1./r_2d**2*drad(r_2d**2*f_visc_r,rr) +\
        1./r_2d/sint_2d*dth(sint_2d*f_visc_t,tt))

savefile = datadir + dirname_stripped + '_torque_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + '.npy'

print('Saving torques at ' + savefile + ' ...')
np.save(savefile, (torque_rs, torque_mc, torque_visc))
