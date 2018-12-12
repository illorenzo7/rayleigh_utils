# Author: Loren Matilsky
# Created on: 05/04/2018
# Computes the radius at the base of the thermal boundary layer by calculating
# where the spherically averaged entropy deviation reaches 10% (can be changed
# with the "tol" variable) its max value

# Import relevant modules here
import numpy as np
from common import get_widest_range_file, interpx
import sys

# Get simulation directory and data directory
dirname = sys.argv[1]
dirname_stripped = dirname.split('/')[-1]
datadir = dirname + '/data/'

s_file = get_widest_range_file(datadir, 's_spherical_mean')
s_spherical_mean = np.load(datadir + s_file)
s_max = np.max(s_spherical_mean)

# Get normalized entropy profile
s_normalized = s_spherical_mean/s_max

# Get basic grid info
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr = len(rr)
# Normalized radius
rrn = rr/ro

tol = 0.1 # Set the bottom of the TBL to where S gets within a fraction tol
    # of its maximum value

ir_close = np.argmin(np.abs(s_normalized - (1 - tol)))
rn_close1 = rrn[ir_close]
rn_close2 = rrn[ir_close + 1]
sn_close1 = s_normalized[ir_close]
sn_close2 = s_normalized[ir_close + 1]
r_tbl = interpx(rn_close1, sn_close1, rn_close2, sn_close2, 1 - tol)

np.save(datadir + dirname_stripped + '_r_tbl.npy', r_tbl)
print('The bottom of the thermal boundary layer is located at %0.3f r_o ... '\
      %r_tbl)