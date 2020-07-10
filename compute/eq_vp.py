# Author: Loren Matilsky
# Created: 03/23/2020
# This script takes the time/longitudinally averaged <v_phi> and outputs
# it in the file "eq_vp" for Rayleigh later to read
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from compute_grid_info import compute_grid_info
from get_parameter import get_parameter
from common import get_widest_range_file, get_dict
dirname = sys.argv[1]

fname = 'eq_vp'
azav_file = None
args = sys.argv[2:]
nargs = len(args)
outdir = dirname
for i in range(nargs):
    arg = args[i]
    if arg == '-usefile':
        azav_file = args[i+1]
    elif arg == '-fname':
        fname = args[i+1]
    elif arg == '-outdir':
        outdir = args[i+1]

# Read in the AZ_Avgs data
datadir = dirname + '/data/'
if azav_file is None:
    azav_file = get_widest_range_file(datadir, 'AZ_Avgs')
print("getting data from data/%s" %azav_file)
di = get_dict(datadir + azav_file)
vals = di['vals']
nt, nr = di['nt'], di['nr']
lut = di['lut']
mean_vp = vals[:, :, lut[3]]

# Write the data
print ("Writing <v_phi> to", outdir + '/' + fname)

f = open(outdir + '/' + fname, 'wb')
sigpi = np.array(314, dtype=np.int32)

f.write(sigpi.tobytes())
for it in range(nt):
    for ir in range(nr):
        f.write(mean_vp[it,ir].tobytes())
f.close()
