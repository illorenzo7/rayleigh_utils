import numpy as np
import matplotlib.pyplot as plt
import sys, os
from diagnostic_reading import ReferenceState
from common import get_widest_range_file, strip_dirname
from get_parameter import get_parameter

dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)
datadir = dirname + '/data/'
plotdir = dirname + '/plots/amom_flux/rslice/'

rs_file = get_widest_range_file(datadir, 'rs')
print ('Reading in Reynolds stress from ' + datadir + rs_file + ' ...')
(vr2_p, vt2_p, vp2_p, vrvp_p, vrvt_p, vtvp_p,\
 vr2_m, vt2_m, vp2_m, vrvp_m, vrvt_m, vtvp_m, fplus, fminus) =\
     np.load(datadir + rs_file)

ref = ReferenceState('reference', dirname + '/')
Hrho = -1/ref.dlnrho

# Get grid info
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr, nt = len(rr), len(tt)
lats = 90 - 180*tt/np.pi 
it1, it2 = np.argmin(np.abs(lats + 15)), np.argmin(np.abs(lats - 15))

# Get frame rotation rate
Om_0 = get_parameter(dirname, 'angular_velocity')

# average velocity amplitudes between +/- 15 degrees.
v2_p = np.mean(((vr2_p + vt2_p + vp2_p)/fplus)[it1:it2, :], axis=0)
v2_m = np.mean(((vr2_m + vt2_m + vp2_m)/fminus)[it1:it2, :], axis=0)

ro_p = np.sqrt(v2_p)/2/Om_0
ro_m = np.sqrt(v2_m)/2/Om_0

plt.plot(rr, ro_p/1e8, label='upflow')
plt.plot(rr, ro_m/1e8, label='downflow')
plt.legend()
