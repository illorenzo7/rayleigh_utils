import numpy as np
import sys

# Get simulation directory and data directory
dirname = sys.argv[1]
dirname_stripped = dirname.split('/')[-1]
datadir = dirname + '/data/'

ir_tbl = np.load(datadir + dirname_stripped + '_ir_tbl.npy')

# Get basic grid info
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr = len(rr)

r_tbl = rr[ir_tbl]
d_tbl_percent = rr_depth[ir_tbl]*100
d_tbl = d_tbl_percent/100*d
print('ir_tbl = %i' %ir_tbl)
print('r_tbl = %.3fro = %.1f Mm' %(r_tbl/ro, r_tbl/1.0e8))
print('depth_tbl = %.1f Mm' %(d_tbl/1.0e8))
print('depth_tbl = %.2f %%' %(d_tbl_percent))