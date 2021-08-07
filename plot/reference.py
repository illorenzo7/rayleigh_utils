# Created: 05/03/2019
# Author: Loren Matilsky

import matplotlib as mpl
import numpy as np
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt

import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])

from common import *
from plotref import plotref

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'

args = sys.argv[2:]
nargs = len(args)
fname = 'equation_coefficients'
xminmax = None
ylog = False
for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-log':
        ylog = True
    elif arg == '-fname':
        fname = args[i+1]
    elif arg == '-crb':
        fname = 'custom_reference_binary'
    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])

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


plt.tight_layout() 
    
axs[0,0].set_title(dirname_stripped, ha='left')
    
plt.savefig(plotdir + dirname_stripped + '_reference_state.png', dpi=300)
plt.show()
