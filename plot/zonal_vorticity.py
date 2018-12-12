####################################################################################################
#
#   Plots the various terms in the gyroscopic pumping equation, performing
#   finite-differences BEFORE taking the time average. I hope this will 
#   yield a smoother balance!
#
#   Loren Matilsky, 11/26/2017
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

import time
import numpy as np
import matplotlib.pyplot as plt
from get_sim_time import *
from derivs import drad, dth
import numpy as np
from azavg_util import *

print('Reading command line arguments... '),
t1 = time.time()
from timeaverage import *
from print_iter_info import *
t2 = time.time()
print('%.2f sec' %(t2 - t1))

t1 = time.time()
data,start_time,end_time,iter_start,iter_end,az=\
        timeaverage(AzAverage)

# Find the times we average between, measured in days (or TDTs) since
# the beginning of the simulation
if (not days):
    t1_since_start = (start_time - simulation_start_time)\
            /thermal_diffusion_time
    t2_since_start = (end_time - simulation_start_time)\
            /thermal_diffusion_time
else:
    t1_since_start = (start_time - simulation_start_time)/86400.
    t2_since_start = (end_time - simulation_start_time)/86400.

# Label the plot using the center of the averaging interval and its length
tcenter = (t1_since_start + t2_since_start)/2.
delta_T = (t2_since_start - t1_since_start)/2.

t2 = time.time()

print('Read data in %.2f sec' %(t2 - t1))


print('Preparing plots... '),
t1 = time.time()
n_r = az.nr
n_t = az.ntheta

# time we averaged over
time_label = r'$\Delta t = %.1f \ \rm{days}$' %delta_T

# variables
##########
# velocities
vr = data[:,:,az.lut[1]]
vt = data[:,:,az.lut[2]]
vphi = data[:,:,az.lut[3]]
vp = vphi

# TD variables
entropy = data[:,:,az.lut[64]]
entropy_dtr = data[:,:,az.lut[88]]

rhovr = data[:,:,az.lut[61]]
rho = rhovr/vr

realnu = nu*2.e12 # 'read_command_line' calculates nu in terms of nu0=2.e12
        # cm^2/sec in order to label the plot
mu = realnu*rho # dynamic viscosity


# radius
radius = az.radius
dradius = np.zeros_like(radius)
dradius[1:n_r] = radius[1:n_r] - radius[:n_r - 1]
dradius[0] = dradius[1]

# angle
sintheta = az.sintheta
costheta = az.costheta
cottheta = costheta/sintheta
theta = np.arccos(costheta)
dtheta = np.zeros_like(theta)
dtheta[1:n_t] = theta[1:n_t] - theta[:n_t - 1]
dtheta[0] = dtheta[1]

vphi_dr = drad(vphi, radius)

vphi_dt = dth(vphi, theta)

#vphi_dr = np.zeros_like(vphi)
#vphi_dr[:,1:n_r] = (vphi[:,1:n_r] - vphi[:,:n_r - 1])/dradius[1:n_r]
#vphi_dr[:,0] = vphi_dr[:,1]

#vphi_dt = np.zeros_like(vphi)
#temp = np.transpose((vphi[1:n_t,:] - vphi[:n_t - 1,:]))/dtheta[1:n_t]
#vphi_dt[1:n_t,:] = np.transpose(temp) 
#vphi_dt[0,:] = vphi_dt[1,:]


vphi_dz = np.zeros_like(vphi)
for it in range(n_t):
    for ir in range(n_r):
        vphi_dz[it,ir] = costheta[it]*vphi_dr[it,ir] - \
                sintheta[it]*vphi_dt[it,ir]/radius[ir]

# terms in thermal wind balance (baroclinic "bc" and stretching term)
# (These are the terms remaining after ignoring viscosity and advection)
newton_G = 6.67259e-8 # cm^3 g^(-1) s^(-2)
local_g = newton_G*central_mass/radius**2

mfb_stretch = np.zeros_like(vphi)
mfb_bc = np.zeros_like(vphi)
for ir in range(n_r):
    mfb_stretch[:,ir] = 2.*angular_velocity*vphi_dz[:,ir]
    mfb_bc[:,ir] = -local_g[ir]/cP*entropy_dtr[:,ir]




# Compute viscous force
#########################
sint_2d = np.zeros((n_t,1))
sint_2d[:,0] = sintheta
cott_2d = np.zeros((n_t,1))
cott_2d[:,0] = cottheta
r_2d = np.zeros((1,n_r))
r_2d[0,:] = radius

# First, rates of strain
s_rr = drad(vr,radius)

s_rt = 0.5*(r_2d*drad(vt/r_2d, radius) + 1./r_2d*dth(vr,theta))

s_tt = 1./r_2d*dth(vt,theta) + 1./r_2d*vr

s_pp = 1./r_2d*vr + 1./r_2d*vt*cott_2d

div_v = s_rr + s_tt + s_pp

# Viscous stress
d_rr = 2.*mu*(s_rr - 1./3.*div_v)
d_rt = 2.*mu*s_rt
d_tt = 2.*mu*(s_tt - 1./3.*div_v)
d_pp = 2.*mu*(s_pp - 1./3.*div_v)

# viscous force
force_r = 1./r_2d**2*drad(r_2d**2*d_rr,radius) +\
        1./(r_2d*sint_2d)*dth(sint_2d*d_rt,theta) -\
        1./r_2d*(d_tt+d_pp)
force_t = 1./r_2d**2*drad(r_2d**2*d_rt,radius) +\
        1./(r_2d*sint_2d)*dth(sint_2d*d_tt,theta) -\
        1./r_2d*cott_2d*d_pp

# viscous term in Taylor_Proudman: (del x visc force)_phi
mfb_visc = 1./r_2d*drad(r_2d/rho*force_t,radius) - \
        1./r_2d*dth(force_r/rho,theta)


# advection/R.S. terms:
udu_r_mean = data[:,:,az.lut[155]]
udu_t_mean = data[:,:,az.lut[156]]


udu_r_fluct = data[:,:,az.lut[164]]

udu_t_fluct = data[:,:,az.lut[165]]


mfb_mc = 1./r_2d*drad(r_2d*udu_t_mean,radius) - \
        1./r_2d*dth(udu_r_mean,theta)
mfb_mc *= -1./rho

mfb_rs = 1./r_2d*drad(r_2d*udu_t_fluct,radius) - \
        1./r_2d*dth(udu_r_fluct,theta)
mfb_rs *= -1./rho

# total "advection term"
mfb_adv = mfb_mc + mfb_rs

# magnitudes
mag_stretch = np.median(np.abs(mfb_stretch))
mag_bc = np.median(np.abs(mfb_bc))
mag_rs = np.median(np.abs(mfb_rs))
mag_adv = np.median(np.abs(mfb_adv))
mag_visc = np.median(np.abs(mfb_visc))
mag_mc = np.median(mfb_mc)

# stds of various terms
stretch_sigma = np.std(mfb_stretch)
bc_sigma = np.std(mfb_bc)
adv_sigma = np.std(mfb_adv)
visc_sigma = np.std(mfb_visc)

rs_sigma = np.std(mfb_rs)
mc_sigma = np.std(mfb_mc)

# maxima of various terms
stretch_maxabs = np.max(np.abs(mfb_stretch))
bc_maxabs = np.max(np.abs(mfb_bc))
adv_maxabs = np.max(np.abs(mfb_adv))
visc_maxabs = np.max(np.abs(mfb_visc))

rs_maxabs = np.max(np.abs(mfb_rs))
mc_maxabs = np.max(np.abs(mfb_mc))


maxabs = nsigma*np.max((stretch_sigma,bc_sigma,visc_sigma,adv_sigma))

my_min = -maxabs
my_max = maxabs

if (user_specified_minmax):
    my_min = user_min
    my_max = user_max

eq_factor = (mfb_stretch + mfb_bc + mfb_adv + mfb_visc)/(np.abs(mfb_stretch) + np.abs(mfb_bc) + np.abs(mfb_adv) + np.abs(mfb_visc))

mag_eq = np.median(np.abs(eq_factor))
eq_sigma = np.std(eq_factor)
eq_min = -nsigma*eq_sigma
eq_max = nsigma*eq_sigma

eq_maxabs = np.max(np.abs(eq_factor))

overallmax = np.max((np.abs(my_min),np.abs(my_max)))
maxabs_exp = np.floor(np.log10(overallmax))
print('maxabs_exp', maxabs_exp)
print('overallmax', overallmax)

mfb_stretch /= 10**maxabs_exp
mfb_bc /= 10**maxabs_exp
mfb_visc /= 10**maxabs_exp
mfb_adv /= 10**maxabs_exp

my_min /= 10**maxabs_exp
my_max /= 10**maxabs_exp


# try to automatically make meaningful ticklabels
my_ticks = np.array([my_min,0.5*my_min, 0., 0.5*my_max, my_max])
my_ticklabels = np.zeros_like(my_ticks, dtype = '|S5')
for ii in range(len(my_ticks)):
    my_ticklabels[ii] = '%1.1f' %my_ticks[ii]

eq_ticks = np.array([eq_min, .5*eq_min, 0., eq_max/2.,eq_max])
eq_ticklabels = np.zeros_like(eq_ticks, dtype = '|S5')
for ii in range(len(eq_ticks)):
    my_ticklabels[ii] = '%.1f' %eq_ticks[ii]

##############################################
#   Here we set up the actual figure.
#   We do a top-row of 3 figures, bottom row of 2 figures
#   Spacing is set manually here

border = 0.1
hspace = 0.1
vspace = 0.1
height = (1.0 - 2*border - vspace)/2
width = (1.0 - 2*border - 2*hspace)/3

if (saveplot):
    fig = plt.figure(figsize=(4.0*3, 10.0*2), dpi=1000)
    plt.rcParams.update({'font.size': 12})
else:
    fig = plt.figure(figsize=(4.0*3, 10.0*2), dpi=80)
    plt.rcParams.update({'font.size': 14})


# Stretching term
tsize = 18
lsize = 14
#ax1 = fig.add_subplot(1,5,1)
units = r'$\ (10^{%i} \ \rm{s}^{-2})$' %maxabs_exp

ax1 =  plt.axes([border, 1.0 - border - height, width, height])
plot_azav(fig,ax1,mfb_stretch,radius,costheta,sintheta,mycmap='RdYlBu_r',\
        boundsfactor = 2., boundstype='manual', units= units,\
        fontsize = lsize,ticks=my_ticks,ticklabels=my_ticklabels, user_minmax=(my_min,my_max),nlevs=my_nlevs,contours=False)
plt.title(r'$2\Omega_0\langle\partial v_\phi/\partial z\rangle$',\
        fontsize=tsize)

ax2 =  plt.axes([border + width + hspace, 1.0 - border - height, width, height])
plot_azav(fig,ax2,mfb_bc,radius,costheta,sintheta,mycmap='RdYlBu_r',\
         boundstype='manual', units=units,\
        fontsize = lsize,ticks=my_ticks,ticklabels=my_ticklabels, user_minmax=(my_min,my_max),nlevs=my_nlevs, contours=False)
plt.title(r'$(g/rc_P)\langle\partial S/\partial\theta\rangle$',\
        fontsize=tsize)

ax3 =  plt.axes([1.0 - border - width, 1.0 - border - height, width, height])

plot_azav(fig,ax3,mfb_adv,radius,costheta,sintheta,mycmap='RdYlBu_r',\
        boundsfactor = 2., boundstype='manual', units=units,\
        fontsize = lsize,ticks=my_ticks,ticklabels=my_ticklabels, user_minmax=(my_min,my_max),nlevs=my_nlevs, contours=False)
plt.title(r'$-\langle\nabla\times (\mathbf{v}\cdot\nabla \mathbf{v})\rangle_\phi$',fontsize = tsize)

ax4 =  plt.axes([border + (width + hspace)/2, border, width, height])
plot_azav(fig,ax4,mfb_visc,radius,costheta,sintheta,mycmap='RdYlBu_r',\
        boundsfactor = 2., boundstype='manual', units=units,\
        fontsize = lsize,ticks=my_ticks,ticklabels=my_ticklabels, user_minmax=(my_min,my_max),nlevs=my_nlevs, contours=False)

plt.title(r'$\langle \nabla\times (\mathbf{f}_{\rm{visc}}/\overline{\rho})\rangle_\phi$',fontsize = tsize)

ax5 =  plt.axes([1.0 - border - width - hspace, border, width, height])
#ax5 =  plt.axes([0.6, 0.0, 0.25, 0.4])
plot_azav(fig,ax5,eq_factor,radius,costheta,sintheta,mycmap='RdYlBu_r',\
        boundsfactor = 2., boundstype='manual', units='',\
        fontsize = lsize,ticks=eq_ticks,ticklabels=eq_ticklabels, user_minmax=(eq_min,eq_max),nlevs=my_nlevs,contours=False)
plt.title(r'$\frac{\rm{total\ vorticity\ flux}}{\rm{vorticity\ flux\ mag.}}$',fontsize = 14)


t2 = time.time()
print('%.2f sec' %(t2 - t1))

#print_iter_info(start_time,end_time)
#plt.subplots_adjust(top=.85)
if (labelplots):
    fig.text(0.5,0.97,generic_title,fontsize=24, \
        horizontalalignment='center',bbox={'facecolor':'white'})
    fig.text(0.5, 0.03,r'$t_c = %.0f, \Delta t = %.1f\ \rm{(days)}$' %(tcenter,delta_T),fontsize=24,\
        horizontalalignment='center',bbox={'facecolor':'white'})

print('median(abs(mfb_stretch)) == %e' %mag_stretch)
print('median(abs(mfb_bc)) == %e' %mag_bc)
print('median(abs(mfb_adv)) == %e' %mag_adv)
print('median(abs(mfb_visc)) == %e' %mag_visc)
print('----------------------------------')
print('median(abs(mfb_rs)) == %e' %mag_rs)
print('median(abs(mfb_mc)) == %e' %mag_mc)
print('-----------------------------------')
print('median(abs(eq_factor)) == %e' %mag_eq)
print('-----------------------------------')
print('maxabs (nisgma*max(sigmas) == %e' %maxabs)
print('max(mfb_stretch) == %e' %stretch_maxabs)
print('max(mfb_bc) == %e' %bc_maxabs)
print('max(mfb_adv) == %e' %adv_maxabs)
print('max(mfb_visc) == %e' %visc_maxabs)
print('max(eq_factor) == %e' %eq_maxabs)
print('----------------------------------')
print('max(mfb_mc) == %e' %mc_maxabs)
print('max(mfb_rs) == %e' %rs_maxabs)

#print('overall max == %e', %(np.max(
print('---------------------------------')
print('nsigma == %.1f' %nsigma)
print('adv_sigma == %e' %adv_sigma)
print('bc_sigma == %e' %bc_sigma)
print('stretch_sigma == %e' %stretch_sigma)
print('visc_sigma == %e' %visc_sigma)
print('eq_sigma == %e' %eq_sigma)




if (saveplot):
    plt.savefig(savefile)  
else:
    plt.show()

