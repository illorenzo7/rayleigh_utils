import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from azavg_util import plot_azav
from binormalized_cbar import MidpointNormalize
from diagnostic_reading import ReferenceState

dirname = sys.argv[1]

datadir = dirname + '/data/'
plotdir = dirname + '/plots/'

if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

ref = ReferenceState(dirname + '/reference', '')
H_rho = -1./ref.dlnrho

# Get grid info
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr, nt = len(rr), len(tt)

H_rho_2d = H_rho.reshape((1, nr)) 

vr2_p,vt2_p,vp2_p,vrvp_p,vrvt_p,vtvp_p,\
    vr2_m,vt2_m,vp2_m, vrvp_m, vrvt_m, vtvp_m, fplus, fminus\
    = np.load(datadir + 'rs_raw.npy') 

vrvp_t = vrvp_m + vrvp_p
vrvt_t = vrvt_m + vrvt_p
vtvp_t = vtvp_m + vtvp_p

vr2_t = vr2_m + vr2_p
vt2_t = vt2_m + vt2_p
vp2_t = vp2_m + vp2_p

# Total velocity
v2_p = vr2_p + vt2_p + vp2_p
v2_m = vr2_m + vt2_p + vp2_m
v2_t = vr2_t + vt2_p + vp2_t

Om = 7.8e-6
ro_p = np.sqrt(v2_p)/(2*Om*H_rho_2d)
ro_m = np.sqrt(v2_m)/(2*Om*H_rho_2d)
ro_t = np.sqrt(v2_t)/(2*Om*H_rho_2d)

ro_p_rslice = np.mean(ro_p,axis=0)
ro_m_rslice = np.mean(ro_m,axis=0)
ro_t_rslice = np.mean(ro_t,axis=0)

plt.plot(rr/ro, ro_p_rslice, label =
            r'upfl3$', color='b')
ax.plot(radius4/rsun, rossby4, label =
            r'$N_\rho=4$', color='orange')
ax.plot(radius5/rsun, rossby5, label =
            r'$N_\rho=5$',color='g')

plt.plot(radius3/rsun, np.ones_like(radius3), 'k--')
# mark the location of the Rossby no. transition for case N4
r_trans4 = 0.9329211
r_trans5 = 0.92187
#plt.yscale('log')
maxes = [np.max(rossby3), np.max(rossby4), np.max(rossby5)]
mins = [np.min(rossby3), np.min(rossby4), np.min(rossby5)]

max_diff = max(maxes) - min(mins)
plt.xlabel(r'$r/R_\odot$',fontsize=14)
plt.ylabel(r'$\rm{Ro_c}\ (v^\prime/2\Omega_0 H_\rho)$',fontsize=14)
plt.xlim((np.min(radius3/rsun),np.max(radius3/rsun)))
plt.ylim((min(mins) - 0.2*max_diff,max(maxes)+0.2*max_diff))
#mark the location of the Ro no. transition
rossby_range = np.linspace(min(mins) - 0.2*max_diff,max(maxes) + 0.2*max_diff,100)
plt.plot(np.ones(100)*r_trans4,rossby_range,linestyle='--',color='orange')
plt.plot(np.ones(100)*r_trans5,rossby_range,'g--')

#plt.title(r'$\rm{Averaged \ over \ %.1f \ days}$'  %delta_T,fontsize=24)

plt.legend(loc='upper left', title='Number of Scale Heights',shadow=True,fontsize=12)    
plt.tight_layout()


plt.plot(rr/ro, ro_p_rslice)
plt.plot(rr/ro, ro_m_rslice)
plt.plot(rr/ro, ro_t_rslice)

plt.show()
