# Created: 05/03/2019
# Author: Loren Matilsky

import matplotlib as mpl
import numpy as np
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}

import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['idref'])

from common import strip_dirname
from plotref import plotref
from rayleigh_diagnostics import ReferenceState, GridInfo
from reference_tools import equation_coefficients

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

try:
    ref = ReferenceState(dirname + '/reference', '')
    r = ref.radius
    T = ref.temperature
    rho = ref.density
    dlnT = ref.dlnt
    dlnrho = ref.dlnrho
    s = ref.entropy
    dsdr = ref.dsdr
    d2lnrho = ref.d2lnrho
    gravity = ref.gravity
    heating = ref.heating
    Q = heating*rho*T
except:
    eq = equation_coefficients()
    eq.read(dirname + '/equation_coefficients')
    r = eq.radius
    T = eq.functions[3]
    rho = eq.functions[0]
    cp = 3.5e8
    gam = 5.0/3.0
    gas_R = (gam - 1.0)*cp/gam
    p = rho*gas_R*T
    dlnT = eq.functions[9]
    dlnrho = eq.functions[7]
    dsdr = eq.functions[13]
    nr = len(dsdr)
    s = np.zeros(nr)
    gi = GridInfo(dirname + '/grid_info')
    rw = gi.rweights
    ri, ro = np.min(r), np.max(r)
    factor = 1.0/3.0*ro**3 - 1.0/3.0*ri**3
    # Remember r is reversed
    r_rev = np.copy(r[::-1])
    rw_rev = np.copy(rw[::-1])
    dsdr_rev = np.copy(dsdr[::-1])
    for ir in range(nr):
        s[nr - 1 - ir] = factor*np.sum((dsdr_rev/r_rev**2*rw_rev)[:ir])     
    d2lnrho = eq.functions[8]

    # Gravity
    buoy = eq.functions[1]
    gravity = buoy*cp/rho

    # Heating
    Q = eq.constants[9]*eq.functions[5]
    
fig, axs = plotref(r, T, rho, p, dlnT, dlnrho, s, dsdr,\
    d2lnrho, gravity, Q, color='k', xminmax=xminmax)

plt.tight_layout() 
    
axs[0,0].set_title('          ' + dirname_stripped, **csfont)
    
plt.savefig(plotdir + dirname_stripped + '_reference_state.png', dpi=300)
plt.show()
