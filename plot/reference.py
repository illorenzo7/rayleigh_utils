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
from get_eq import get_eq

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Get other arguments
xminmax = None
rvals = []
ylog = True

args = sys.argv[2:]
nargs = len(args)
fname = 'equation_coefficients'
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
        fname = args[i+1]
    elif arg == '-crb':
        fname = 'custom_reference_binary'
    elif arg == '-nolog':
        ylog = False

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

eq = get_eq(dirname, fname=fname)
r = eq.radius
T = eq.temperature
rho = eq.density
p = eq.pressure
dlnT = eq.dlnT
dlnrho = eq.dlnrho
dsdr = eq.dsdr
# Integrate to obtain s(r)
s = cumtrapz(dsdr, r, initial=0)
d2lnrho = eq.d2lnrho
gravity = eq.gravity
Q = eq.Q
    
fig, axs = plotref(r, T, rho, p, dlnT, dlnrho, s, dsdr,\
    d2lnrho, gravity, Q, color='k', xminmax=xminmax, ylog=ylog)

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
    
axs[0,0].set_title(dirname_stripped, ha='left', **csfont)
    
plt.savefig(plotdir + dirname_stripped + '_reference_state.png', dpi=300)
plt.show()
