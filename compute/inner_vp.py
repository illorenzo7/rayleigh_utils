# Author: Loren Matilsky
# Created: 03/23/2020
# This script takes the time/longitudinally averaged <v_phi> at the 
# base of the CZ for an [equilibrated] RZ-CZ run and outputs <v_phi>
# vs. theta in the file "inner_vp" for Rayleigh later to read
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from compute_grid_info import compute_grid_info
from get_parameter import get_parameter
from common import get_widest_range_file, get_dict
dirname = sys.argv[1]

fname = 'inner_vp'
azav_file = None
args = sys.argv[2:]
nargs = len(args)
two_domains = True
r0 = None
for i in range(nargs):
    arg = args[i]
    if arg == '-usefile':
        azav_file = args[i+1]
    elif arg == '-cz':
        two_domains = False
    elif arg == '-fname':
        fname = args[i+1]
    elif arg == '-r0':
        r0 = float(args[i+1])

# Read in the AZ_Avgs data
datadir = dirname + '/data/'
if azav_file is None:
    azav_file = get_widest_range_file(datadir, 'AZ_Avgs')
print("getting data from data/", azav_file)
di = get_dict(datadir + azav_file)
vals = di['vals']
lut = di['lut']
rr = di['rr']
mean_vp = vals[:, :, lut[3]]

# Get relevant info from main_input file
if two_domains:
    domain_bounds = tuple(get_parameter(dirname, 'domain_bounds'))
    ncheby = tuple(get_parameter(dirname, 'ncheby'))
    ir_base_CZ = ncheby[1] - 1
else:
    nr = get_parameter(dirname, 'n_r')
    ir_base_CZ = nr - 1

if r0 is None:
    ir0 = ir_base_CZ
else:
    ir0 = np.argmin(np.abs(rr - r0))
    
inner_vp = mean_vp[:, ir0]

# Write the data
print ("Writing <v_phi> as a function of theta for ")
print("r0 = %1.3e cm" %rr[ir0])

f = open(dirname + '/' + fname, 'wb')
sigpi = np.array(314, dtype=np.int32)

f.write(sigpi.tobytes())
f.write(inner_vp.tobytes())
f.close()
