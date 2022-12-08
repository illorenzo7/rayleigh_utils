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
kwargs_default = dotdict({'radtype': 'sslice'})
kw = update_dict(kwargs_default, clas)
radtype = kw.radtype

if radtype == 'sslice':
    dataname = 'Shell_Slices'
    samplelabel = 'rvals'
if radtype == 'merslice':
    dataname = 'Meridional_Slices'
    samplelabel = 'lonvals'
if radtype == 'eqslice':
    dataname = 'Equatorial_Slices'

# state what we are about to print
print (buff_line)
print (make_bold("for     %20s:" %dirname_stripped))
print (make_bold("for the %20s:" %dataname))
# get the sampled values
di = get_sliceinfo(dirname, dataname=dataname)
print (buff_line)
print ("%i variables sampled" %di.nsamplevals)
print ('qv =', arr_to_str(di.qv, '%i'))

if not radtype == 'eqslice':
    print (buff_line)
    print ("%i sampling locations" %di.nsamplevals)
    print ('iisamplevals     =', arr_to_str(np.arange(di.nsamplevals), '%9i'))
    print ('isamplevals      =', arr_to_str(di.isamplevals, '%9i'))
    print (buff_line)
    print ('samplevals       = ' + arr_to_str(di.samplevals, "%9.3e"))
    print (buff_line)
