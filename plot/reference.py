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
