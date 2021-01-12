# NEEDS UPDATED DOCS AND TESTING
#   Author: Loren Matilsky
#   Created: 11/26/2017
#  
#   Plots the various terms in the gyroscopic pumping equation, performing
#   finite-differences BEFORE taking the time average. I hope this will 
#   yield a smoother balance!
#
#
#   This plotting routine makes use of the azimuthal averages class, which 
#   has the following attributes:
#    ----------------------------------
#    self.niter                                    : number of time steps
#    self.nq                                       : number of diagnostic quantities output
#    self.nr                                       : number of radial points
#    self.ntheta                                   : number of theta points
#    self.qv[0:nq-1]                               : quantity codes for the diagnostics output
#    self.radius[0:nr-1]                           : radial grid
#    self.costheta[0:ntheta-1]                     : cos(theta grid)
#    self.sintheta[0:ntheta-1]                     : sin(theta grid)
#    self.vals[0:ntheta-1,0:nr-1,0:nq-1,0:niter-1] : The phi-averaged diagnostics 
#    self.iters[0:niter-1]                         : The time step numbers stored in this output file
#    self.time[0:niter-1]                          : The simulation time corresponding to each time step
#    self.version                                  : The version code for this particular output (internal use)
#    self.lut                                      : Lookup table for the different diagnostics output
#    ----------------------------------

import numpy as np
import os, sys
from derivs import drad, dth
import numpy as np
from diagnostic_reading import AzAverage, ReferenceState
from get_parameter import get_parameter

dirname = sys.argv[1]
radatadir = dirname + '/AZ_Avgs/'
ref = ReferenceState(dirname + '/reference','')
angular_velocity = get_parameter(dirname, 'angular_velocity')
nu_top = get_parameter(dirname, 'nu_top')
nu = nu_top
central_mass = get_parameter(dirname, 'poly_mass')
cP = get_parameter(dirname, 'pressure_specific_heat')

datadir = dirname + '/data/'
if (not os.path.isdir(datadir)):
    os.makedirs(datadir)

files = os.listdir(radatadir)
nfiles = len(files)
files.sort()

times,dtimes,weights = np.load(datadir + 'azavgs_times.npy')

az0 = AzAverage(radatadir + files[0], '')
nr = az0.nr
nt = az0.ntheta
rr = az0.radius
cost = az0.costheta
sint = az0.sintheta
rr_2d = rr.reshape((1,nr))

tt = np.arccos(cost)
tt_2d = tt.reshape((nt,1))
cost_2d = np.cos(tt_2d)
sint_2d = np.sin(tt_2d)
cott_2d = cost_2d/sint_2d

# grid spacings

# radius
dr = np.zeros_like(rr)
dr[1:nr] = rr[1:nr] - rr[:nr - 1]
dr[0] = dr[1]

# angle
dt = np.zeros_like(tt)
dt[1:nt] = tt[1:nt] - tt[:nt - 1]
dt[0] = dt[1]

vavg = np.load(datadir + 'vavg.npy')
vr_av, vt_av, vp_av = vavg

rho = ref.density
rho_2d = rho.reshape((1,nr))
mu = nu*rho_2d # dynamic viscosity

gp_stretch = np.zeros_like(vp_av)
gp_visc = np.zeros_like(vp_av)
gp_bc = np.zeros_like(vp_av)

gp_rs = np.zeros_like(vp_av)
gp_mc = np.zeros_like(vp_av)
gp_adv = np.zeros_like(vp_av)

count = 0
#for ii in range(nfiles):
for ii in range(100): # for debugging purposes
    az = AzAverage(radatadir + files[ii], '')
    print ('adding  azavg %s to the average...' %files[ii])
    niter = az.niter
    for tindex in range(niter):
        # velocities
        vr = az.vals[:,:,az.lut[1],tindex]
        vt = az.vals[:,:,az.lut[2],tindex]
        vp = az.vals[:,:,az.lut[3],tindex]

        # TD variables
        entropy = az.vals[:,:,az.lut[64],tindex]
        entropy_dtr = az.vals[:,:,az.lut[88],tindex]



        vp_dr = drad(vp, rr)
        vp_dt = dth(vp, tt)

        vp_dz = np.zeros_like(vp_av)
        for it in range(nt):
            for ir in range(nr):
                vp_dz[it,ir] = cost[it]*vp_dr[it,ir] - \
                        sint[it]*vp_dt[it,ir]/rr[ir]

        # terms in thermal wind balance (baroclinic "bc" and stretching term)
        # (These are the terms remaining after ignoring viscosity and advection)
        newton_G = 6.67259e-8 # cm^3 g^(-1) s^(-2)
        local_g = newton_G*central_mass/rr**2

        mfb_stretch = np.zeros_like(vp)
        mfb_bc = np.zeros_like(vp)
        for ir in range(nr):
            mfb_stretch[:,ir] = 2.*angular_velocity*vp_dz[:,ir]
            mfb_bc[:,ir] = -local_g[ir]/cP*entropy_dtr[:,ir]




        # Compute viscous force
        #########################

        # First, rates of strain
        s_rr = drad(vr,rr)

        s_rt = 0.5*(rr_2d*drad(vt/rr_2d, rr) + 1./rr_2d*dth(vr,tt))

        s_tt = 1./rr_2d*dth(vt,tt) + 1./rr_2d*vr

        s_pp = 1./rr_2d*vr + 1./rr_2d*vt*cott_2d

        div_v = s_rr + s_tt + s_pp

        # Viscous stress
        d_rr = 2.*mu*(s_rr - 1./3.*div_v)
        d_rt = 2.*mu*s_rt
        d_tt = 2.*mu*(s_tt - 1./3.*div_v)
        d_pp = 2.*mu*(s_pp - 1./3.*div_v)

        # viscous force
        force_r = 1./rr_2d**2*drad(rr_2d**2*d_rr,rr) +\
                1./(rr_2d*sint_2d)*dth(sint_2d*d_rt,tt) -\
                1./rr_2d*(d_tt+d_pp)
        force_t = 1./rr_2d**2*drad(rr_2d**2*d_rt,rr) +\
                1./(rr_2d*sint_2d)*dth(sint_2d*d_tt,tt) -\
                1./rr_2d*cott_2d*d_pp

        # viscous term in Taylor_Proudman: (del x visc force)_phi
        mfb_visc = 1./rr_2d*drad(rr_2d/rho_2d*force_t,rr) - \
                1./rr_2d*dth(force_r/rho_2d,tt)


        # advection/R.S. terms:
        udu_r_mean = az.vals[:,:,az.lut[155],tindex]
        udu_t_mean = az.vals[:,:,az.lut[156],tindex]


        udu_r_fluct = az.vals[:,:,az.lut[164],tindex]

        udu_t_fluct = az.vals[:,:,az.lut[165],tindex]


        mfb_mc = 1./rr_2d*drad(rr_2d*udu_t_mean,rr) - \
                1./rr_2d*dth(udu_r_mean,tt)
        mfb_mc *= -1./rho

        mfb_rs = 1./rr_2d*drad(rr_2d*udu_t_fluct,rr) - \
                1./rr_2d*dth(udu_r_fluct,tt)
        mfb_rs *= -1./rho

        # total "advection term"
        mfb_adv = mfb_mc + mfb_rs

        gp_stretch += mfb_stretch*weights[count]
        gp_bc += mfb_bc*weights[count]
        gp_visc += mfb_visc*weights[count]

        gp_rs += mfb_rs*weights[count]
        gp_mc += mfb_mc*weights[count]
        gp_adv += mfb_adv*weights[count]

        count += 1

savefile = datadir + 'gp_terms.npy'
np.save(savefile, (gp_stretch,gp_bc,gp_visc, gp_rs,gp_mc,gp_adv))
