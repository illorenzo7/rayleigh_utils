import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
from common import *
from cla_util import *

# Get CLAs
args = sys.argv 
clas0, clas = read_clas(args)
dirname = clas0.dirname
dirname_stripped = strip_dirname(dirname)

# SPECIFIC ARGS
kwargs_default = dotdict(dict({'type': 'moll'}))
kw = update_dict(kwargs_default, clas)

if kw.type == 'moll':
    datatype = 'Shell_Slices'
if kw.type == 'speclm':
    datatype = 'Shell_Spectra'

# get the radial levels
radlevs = get_slice_levels(dirname, datatype=datatype)
print ('irvals     =', radlevs.inds)
print ('rvals      = ' + arr_to_str(radlevs.radius, "%1.3e"))
print ('rvals/rsun = ' + arr_to_str(radlevs.radius/rsun, "%.3f"))
print (buff_line)
