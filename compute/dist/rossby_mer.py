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

# Plot radial angular momentum transport

fig, ax = plt.subplots()
plot_azav(fig, ax, ro_m, rr, cost, sint,
        contours=False, notfloat=False, units='')
plt.title(r'$({\rm{Ro}}_{\rm{c}})_+$',fontsize=16)
plt.tight_layout()
plt.show()
plt.savefig(plotdir + 'rossby_mer_p.png')
plt.close()

