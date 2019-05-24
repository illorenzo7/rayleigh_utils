# Created: 05/03/2019
# Author: Loren Matilsky

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
#import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
from get_parameter import get_parameter

from common import strip_dirname
from plotref import plotref
from rayleigh_diagnostics import TransportCoeffs
# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

t = TransportCoeffs(dirname + '/transport', '')
r = t.radius
nu = t.nu
dlnu = t.dlnu
k = t.kappa
dlnk = t.dlnkappa

magnetism = get_parameter(dirname, 'magnetism')
if magnetism:
    eta = t.eta
    dlneta = t.dlneta

fig = plt.figure(
fig, axs = plotref(r, T, rho, p, dlnT, dlnrho, dlnp, s, dsdr,\
    d2lnrho, color='k')

plt.tight_layout() 
    
axs[0,0].set_title('          ' + dirname_stripped, **csfont)
    
plt.savefig(plotdir + dirname_stripped + '_reference_state.png', dpi=300)
plt.show()
