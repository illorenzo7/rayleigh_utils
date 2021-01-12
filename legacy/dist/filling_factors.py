import numpy as np
import matplotlib.pyplot as plt
import sys, os
from diagnostic_reading import ReferenceState
from common import get_widest_range_file, strip_dirname

dirname = sys.argv[1]
# Get the stripped name to use in file naming
dirname_stripped = strip_dirname(dirname)

datadir = dirname + '/data/'
plotdir = dirname + '/plots/'

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

rs_file = get_widest_range_file(datadir, 'rs')
print ('Reading in Reynolds stress from ' + datadir + rs_file + ' ...')
vr2_p, vt2_p, vp2_p, vrvp_p, vrvt_p, vtvp_p, vr2_m, vt2_m, vp2_m, vrvp_m, vrvt_m, vtvp_m, fplus, fminus = np.load(datadir + rs_file)

sint_2d = sint.reshape((nt, 1))
fplus_r = np.sum(fplus*sint_2d, axis=0)/np.sum(sint_2d)
fminus_r = np.sum(fminus*sint_2d, axis=0)/np.sum(sint_2d)

plt.plot(rr/ro, fplus_r, label='up')
plt.plot(rr/ro, fminus_r, label='down')
plt.legend()
plt.xlim((ri/ro, 1))
plt.ylim((0, 1))
plt.xlabel(r'$r/r_o$')
plt.ylabel('filling factor')
plt.title(dirname_stripped, fontsize=18)
# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')


plt.savefig(plotdir + dirname_stripped + '_filling_factors.png', dpi=300)
plt.close()
