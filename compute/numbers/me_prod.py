# Author: Loren Matilsky
# Created: 04/04/2021
# This script computes the volume-averaged mag. production terms for a 
# Rayleigh run in directory [dirname]
# Since B_phi is usually by far the strongest component, this is really
# only useful for examining toroidal field generation
# Displays the production terms (and ratio to total induction) 
# at the terminal
# computes in two zones if user specifies, e.g.
# -rbcz 0.718

import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import Shell_Avgs, GridInfo
from common import *

# Get directory name
dirname = sys.argv[1]

# set defaults and adjust them from command line
datadir = dirname + '/data/'
the_file = get_widest_range_file(datadir, 'Shell_Avgs')
rbcz = None
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-usefile':
        the_file = args[i+1]
        the_file = the_file.split('/')[-1]
    if arg == '-rbcz':
        rbcz = float(args[i+1])

# Read in the Shell_Avgs data
print ('ME production terms from ' + datadir + the_file)
di = get_dict(datadir + the_file)
vals = di['vals']
lut = di['lut']
rr = di['rr']
nr = di['nr']

# Read in grid info for radial weights and reference velocity
gi = GridInfo(dirname + '/grid_info')
rw = gi.rweights

# get mag. energy production terms
indtot = vals[:, lut[2019]]
shear = vals[:, lut[2025]]
adv = vals[:, lut[2026]]
comp = vals[:, lut[2027]]
diff = vals[:, lut[2043]]

if not rbcz is None:
    irbcz = np.argmin(np.abs(rr/rsun - rbcz))
    if not (irbcz == 0 or irbcz == nr - 1):
        rwcz = rw[:irbcz+1]/np.sum(rw[:irbcz+1])
        rwrz = rw[irbcz+1:]/np.sum(rw[irbcz+1:])
    else:
        print ('nonD_numbers(): dude, you entered a stupid value for')
        print ('rbcz. you set rbcz = %1.3e' %rbcz)
        print ('it needs be in the (exclusive) range (%.3f, %.3f)' %(np.min(rr)/rsun, np.max(rr)/rsun))
        print ('resetting rbcz = None')
        rbcz = None

# Make empty dictionary
di_out = dict([])
di_out['indtot'] = indtot
di_out['shear '] = shear
di_out['adv   '] = adv
di_out['comp  '] = comp
di_out['diff  '] = diff

base_keys = list(di_out.keys())
for key in base_keys:
    di_out[key + '_gav'] = np.sum(di_out[key]*rw)
    if not rbcz is None:
        di_out[key + '_cz'] = np.sum(di_out[key][:irbcz+1]*rwcz)
        di_out[key + '_rz'] = np.sum(di_out[key][irbcz+1:]*rwrz)

# print everything
print ("===================================")
for key in base_keys:
    exts = ['gav']
    if not rbcz is None:
        exts.append('cz')
        exts.append('rz')
    for ext in exts:
        #print (key + ' (%4s, cgs):     %1.3e'\
        #        %(ext, di_out[key + '_' + ext]))
        print (key + ' (%4s, /indtot): %1.3e'\
                %(ext, di_out[key + '_' + ext]/di_out['indtot_' + ext]))
    print ("===================================")
