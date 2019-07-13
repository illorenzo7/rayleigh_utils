# Created: 05/03/2019
# Author: Loren Matilsky

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}

import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
sys.path.append(os.environ['idref'])

from common import strip_dirname
from plotref import plotref
from rayleigh_diagnostics import ReferenceState
# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Get other arguments
xminmax = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

ref = ReferenceState(dirname + '/reference', '')
r = ref.radius
T = ref.temperature
rho = ref.density
p = ref.pressure
dlnT = ref.dlnt
dlnrho = ref.dlnrho
dlnp = dlnT + dlnrho
s = ref.entropy
dsdr = ref.dsdr
d2lnrho = ref.d2lnrho

fig, axs = plotref(r, T, rho, p, dlnT, dlnrho, dlnp, s, dsdr,\
    d2lnrho, color='k', xminmax=xminmax)

plt.tight_layout() 
    
axs[0,0].set_title('          ' + dirname_stripped, **csfont)
    
plt.savefig(plotdir + dirname_stripped + '_reference_state.png', dpi=300)
plt.show()
