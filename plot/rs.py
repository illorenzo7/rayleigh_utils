import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
from azavg_util import plot_azav
from binormalized_cbar import MidpointNormalize
from diagnostic_reading import ReferenceState

dirname = sys.argv[1]

datadir = dirname + '/data/'
plotdir = dirname + '/plots/rs/'

if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

ref = ReferenceState(dirname + '/reference', '')
rho = ref.density

# Get grid info
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr, nt = len(rr), len(tt)

rho_2d = rho.reshape((1, nr))
r_2d = rr.reshape((1,nr))
sint_2d = sint.reshape((nt, 1))
rho_r_sint = rho_2d*r_2d*sint_2d

vr2_p,vt2_p,vp2_p,vrvp_p,vrvt_p,vtvp_p,\
    vr2_m,vt2_m,vp2_m, vrvp_m, vrvt_m, vtvp_m, fplus, fminus\
    = np.load(datadir + 'rs_pm.npy') 

vrvp_t = vrvp_m + vrvp_p
vrvt_t = vrvt_m + vrvt_p
vtvp_t = vtvp_m + vtvp_p

vr2_t = vr2_m + vr2_p
vt2_t = vt2_m + vt2_p
vp2_t = vp2_m + vp2_p

# Plot radial angular momentum transport
flux_bound_r = max( np.max(np.abs(rho_r_sint*vrvp_p)), np.max(np.abs(rho_r_sint*vrvp_m)),np.max(np.abs(rho_r_sint*vrvp_t)) )

flux_bound_t = max( np.max(np.abs(rho_r_sint*vtvp_p)), np.max(np.abs(rho_r_sint*vtvp_m)),np.max(np.abs(rho_r_sint*vtvp_t)) )

fig, ax = plt.subplots()
plot_azav(fig, ax, rho_r_sint*vrvp_p, rr, cost, sint, \
        boundstype='manual', caller_minmax=(-flux_bound_r, flux_bound_r),
        plotcontours=False, units=r'$\frac{\rm{g}}{\rm{s}^2}$')
plt.title(r'$\rho r\ \sin\theta\ \langle v_r^\prime v_\phi^\prime\rangle_+$',fontsize=16)
plt.tight_layout()
plt.savefig(plotdir + 'amom_flux_r_p.png')
plt.close()

fig, ax = plt.subplots()
plot_azav(fig, ax, rho_r_sint*vrvp_m, rr, cost, sint,\
        plotcontours=False,boundstype='manual', caller_minmax=(-flux_bound_r, flux_bound_r), units=r'$\frac{\rm{g}}{\rm{s}^2}$')
plt.title(r'$\rho r\ \sin\theta\ \langle v_r^\prime v_\phi^\prime\rangle_-$',fontsize=16)
plt.tight_layout()
plt.savefig(plotdir + 'amom_flux_r_m.png')
plt.close()

fig, ax = plt.subplots()
plot_azav(fig, ax, rho_r_sint*vrvp_t, rr, cost, sint,\
        plotcontours=False,boundstype='manual', caller_minmax=(-flux_bound_r, flux_bound_r), units=r'$\frac{\rm{g}}{\rm{s}^2}$')
plt.title(r'$\rho r\ \sin\theta\ \langle v_r^\prime v_\phi^\prime\rangle_t$',fontsize=16)
plt.tight_layout()
plt.savefig(plotdir + 'amom_flux_r_t.png')
plt.close()

# Plot radial angular momentum transport
fig, ax = plt.subplots()
plot_azav(fig, ax, rho_r_sint*vtvp_p, rr, cost, sint, \
        plotcontours=False,boundstype='manual', caller_minmax=(-flux_bound_t, flux_bound_t), units=r'$\frac{\rm{g}}{\rm{s}^2}$')
plt.title(r'$\rho r\ \sin\theta\ \langle v_\theta^\prime v_\phi^\prime\rangle_+$',fontsize=16)
plt.tight_layout()
plt.savefig(plotdir + 'amom_flux_t_p.png')
plt.close()

fig, ax = plt.subplots()
plot_azav(fig, ax, rho_r_sint*vtvp_m, rr, cost, sint,\
        plotcontours=False,boundstype='manual', caller_minmax=(-flux_bound_t, flux_bound_t), units=r'$\frac{\rm{g}}{\rm{s}^2}$')
plt.title(r'$\rho r\ \sin\theta\ \langle v_\theta^\prime v_\phi^\prime\rangle_-$',fontsize=16)
plt.tight_layout()
plt.savefig(plotdir + 'amom_flux_t_m.png')
plt.close()

fig, ax = plt.subplots()
plot_azav(fig, ax, rho_r_sint*vtvp_t, rr, cost, sint,\
        plotcontours=False,boundstype='manual', caller_minmax=(-flux_bound_t, flux_bound_t), units=r'$\frac{\rm{g}}{\rm{s}^2}$')
plt.title(r'$\rho r\ \sin\theta\ \langle v_\theta^\prime v_\phi^\prime\rangle_t$',fontsize=16)
plt.tight_layout()
plt.savefig(plotdir + 'amom_flux_t_t.png')
plt.close()

# Plot the Reynolds stress correlations

# vrvp
fig, ax = plt.subplots()
cor = vrvp_p/np.sqrt(vr2_p*vp2_p)
std_cor = np.std(cor)
plot_azav(fig, ax, cor, rr, cost, sint, \
        plotcontours=False, boundstype='manual', caller_minmax=(-1,1),  units='')
plt.title(r'$(\langle v_r^\prime v_\phi^\prime\rangle/(\langle v_r^{\prime 2}\rangle \langle v_\phi^{\prime 2}\rangle)^{1/2})_+$',\
                fontsize=16)
plt.tight_layout()
plt.savefig(plotdir + 'cor_vrvp_p.png')
plt.close()

fig, ax = plt.subplots()
cor = vrvp_m/np.sqrt(vr2_m*vp2_m)
std_cor = np.std(cor)
plot_azav(fig, ax, cor, rr, cost, sint, \
        plotcontours=False, boundstype='manual', caller_minmax=(-1,1),  units='')
plt.title(r'$(\langle v_r^\prime v_\phi^\prime\rangle/(\langle v_r^{\prime 2}\rangle \langle v_\phi^{\prime 2}\rangle)^{1/2})_-$',\
                fontsize=16)
plt.tight_layout()
plt.savefig(plotdir + 'cor_vrvp_m.png')
plt.close()

fig, ax = plt.subplots()
cor = vrvp_t/np.sqrt(vr2_t*vp2_t)
std_cor = np.std(cor)
plot_azav(fig, ax, cor, rr, cost, sint, \
        plotcontours=False, boundstype='manual', caller_minmax=(-1,1),  units='')
plt.title(r'$(\langle v_r^\prime v_\phi^\prime\rangle/(\langle v_r^{\prime 2}\rangle \langle v_\phi^{\prime 2}\rangle)^{1/2})_t$',\
                fontsize=16)
plt.tight_layout()
plt.savefig(plotdir + 'cor_vrvp_t.png')
plt.close()

# vtvp
fig, ax = plt.subplots()
cor = vtvp_p/np.sqrt(vt2_p*vp2_p)
std_cor = np.std(cor)
plot_azav(fig, ax, cor, rr, cost, sint, \
        plotcontours=False, boundstype='manual', caller_minmax=(-1,1),  units='')
plt.title(r'$(\langle v_\theta^\prime v_\phi^\prime\rangle/(\langle v\theta^{\prime 2}\rangle \langle v_\phi^{\prime 2}\rangle)^{1/2})_+$',\
                fontsize=16)
plt.tight_layout()
plt.savefig(plotdir + 'cor_vtvp_p.png')
plt.close()

fig, ax = plt.subplots()
cor = vtvp_m/np.sqrt(vt2_m*vp2_m)
std_cor = np.std(cor)
plot_azav(fig, ax, cor, rr, cost, sint, \
        plotcontours=False, boundstype='manual', caller_minmax=(-1,1),  units='')
plt.title(r'$(\langle v_\theta^\prime v_\phi^\prime\rangle/(\langle v_\theta^{\prime 2}\rangle \langle v_\phi^{\prime 2}\rangle)^{1/2})_-$',\
                fontsize=16)
plt.tight_layout()
plt.savefig(plotdir + 'cor_vtvp_m.png')
plt.close()

fig, ax = plt.subplots()
cor = vtvp_t/np.sqrt(vt2_t*vp2_t)
std_cor = np.std(cor)
plot_azav(fig, ax, cor, rr, cost, sint, \
        plotcontours=False, boundstype='manual', caller_minmax=(-1,1),  units='')
plt.title(r'$(\langle v_\theta^\prime v_\phi^\prime\rangle/(\langle v_\theta^{\prime 2}\rangle \langle v_\phi^{\prime 2}\rangle)^{1/2})_t$',\
                fontsize=16)
plt.tight_layout()
plt.savefig(plotdir + 'cor_vtvp_t.png')
plt.close()

# vrvt
fig, ax = plt.subplots()
cor = vrvt_p/np.sqrt(vr2_p*vt2_p)
std_cor = np.std(cor)
plot_azav(fig, ax, cor, rr, cost, sint, \
        plotcontours=False, boundstype='manual', caller_minmax=(-1,1),  units='')
plt.title(r'$(\langle v_r^\prime v_\theta^\prime\rangle/(\langle v_r^{\prime 2}\rangle \langle v_\theta^{\prime 2}\rangle)^{1/2})_+$',\
                fontsize=16)
plt.tight_layout()
plt.savefig(plotdir + 'cor_vrvt_p.png')
plt.close()

fig, ax = plt.subplots()
cor = vrvt_m/np.sqrt(vr2_m*vt2_m)
std_cor = np.std(cor)
plot_azav(fig, ax, cor, rr, cost, sint, \
        plotcontours=False, boundstype='manual', caller_minmax=(-1,1),  units='')
plt.title(r'$(\langle v_r^\prime v_\theta^\prime\rangle/(\langle v_r^{\prime 2}\rangle \langle v_\theta^{\prime 2}\rangle)^{1/2})_-$',\
                fontsize=16)
plt.tight_layout()
plt.savefig(plotdir + 'cor_vrvt_m.png')
plt.close()

fig, ax = plt.subplots()
cor = vrvt_t/np.sqrt(vr2_t*vt2_t)
std_cor = np.std(cor)
plot_azav(fig, ax, cor, rr, cost, sint, \
        plotcontours=False, boundstype='manual', caller_minmax=(-1,1),  units='')
plt.title(r'$(\langle v_r^\prime v_\theta^\prime\rangle/(\langle v_r^{\prime 2}\rangle \langle v_\theta^{\prime 2}\rangle)^{1/2})_t$',\
                fontsize=16)
plt.tight_layout()
plt.savefig(plotdir + 'cor_vrvt_t.png')
plt.close()

