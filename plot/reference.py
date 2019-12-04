# Created: 05/03/2019
# Author: Loren Matilsky

import matplotlib as mpl
import numpy as np
from scipy.integrate import cumtrapz
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}

import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['idref'])

from common import strip_dirname, rsun
from plotref import plotref
from rayleigh_diagnostics import ReferenceState, GridInfo
from reference_tools import equation_coefficients

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Get other arguments
xminmax = None
rvals = None

args = sys.argv[2:]
nargs = len(args)
custom_name = None
for i in range(nargs):
    arg = args[i]
    if arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for rval_str in rvals_str:
            rvals.append(float(rval_str))
    elif arg == '-fname':
        custom_name = args[i+1]
    elif arg == '-crb':
        custom_name = 'custom_reference_binary'

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
    if custom_name is None:
        custom_name = 'equation_coefficients'
    eq.read(dirname + '/' + custom_name)
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

    # Integrate to obtain s(r)
    s = cumtrapz(dsdr, r, initial=0)
    d2lnrho = eq.functions[8]

    # Gravity
    buoy = eq.functions[1]
    gravity = buoy*cp/rho

    # Heating
    Q = eq.constants[9]*eq.functions[5]
    
fig, axs = plotref(r, T, rho, p, dlnT, dlnrho, s, dsdr,\
    d2lnrho, gravity, Q, color='k', xminmax=xminmax)

# Mark radii if desired
if not rvals is None:
    for ax in axs.flatten():
        ymin, ymax = ax.get_ylim()
        yvals = np.linspace(ymin, ymax, 100)
        for rval in rvals:
            rval_n = rval/rsun
    #        plt.ylim(ymin, ymax)
            ax.plot(rval_n + np.zeros(100), yvals, 'k--', linewidth=0.8)

plt.tight_layout() 
    
axs[0,0].set_title('          ' + dirname_stripped, **csfont)
    
plt.savefig(plotdir + dirname_stripped + '_reference_state.png', dpi=300)
plt.show()
