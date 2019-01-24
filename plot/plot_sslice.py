import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['co'])
sys.path.append(os.environ['rapp'])
from binormalized_cbar import MidpointNormalize
from mpl_toolkits.basemap import Basemap, addcyclic
from matplotlib import colors
from varprops import texlabels, texunits, var_indices, var_indices_old
from common import get_widest_range_file, get_file_lists
from rayleigh_diagnostics import Shell_Slices
from get_sslice import get_sslice
from sslice_util import show_ortho

dirn = sys.argv[1]
args = sys.argv[2:]
nargs = len(args)

file_list, int_file_list, nfiles = get_file_lists(dirn + '/Shell_Slices/')
fname = file_list[-1] # By default plot the last shell slice
var = 'vr' # By default plot the radial velocity
idepth = 0 # By default plot the shell slice closest to surface
for i in range(nargs):
    arg = args[i]
    if arg == '-n':
        desired_iter = int(args[i+1])
        fname = file_list[np.argmin(np.abs(int_file_list - desired_iter))]
        print ("Couldn't find Shell_Slice # %s, using %s instead"\
               %(str(desired_iter).zfill(8), fname))
    elif arg == '-var':
        var = args[i+1]
    elif arg == '-idepth':
        idepth = int(args[i+1])
        
a = Shell_Slices(fname, dirn + '/Shell_Slices/')

sslice = get_sslice(a, var, dirname=dirn)

ri = 5.
ro = 6.586
rloc = 6.586
try:
    rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(dirn + '/data/grid_info.npy')
    rloc = a.radius[idepth]
except:
    print ('No grid_info.npy')
    
show_ortho(sslice[:, :, idepth], var, ri=ri, ro=ro, rloc=rloc, iterr=a.iters[0])