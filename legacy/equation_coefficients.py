# Created: 08/05/2019
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

from common import *
from plotref import plotref
from reference_tools import equation_coefficients
from get_parameter import get_parameter
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

eq = equation_coefficients()
eq.read(dirname + '/equation_coefficients')

# radius
r = eq.radius
nr = eq.nr

# thermodynamic variables, including pressure
T = eq.functions[3, :]
rho = eq.functions[0, :]
cp = get_parameter(dirname, 'pressure_specific_heat')
gam = 5/3
R = (gam - 1.)*cp/gam
p = rho*R*T

# heating
heating = eq.functions[5, :]/(rho*T)

# log derivatives
dlnrho = eq.functions[7, :]
d2lnrho = eq.functions[8, :]
dlnT = eq.functions[9, :]
dlnp = dlnT + dlnrho

dsdr = eq.functions[13, :]
# compute s from integral of dsdr:
s = np.zeros(nr) 
for i in range(1, nr):
    s[i] = np.trapz(dsdr[:i+1], x=r[:i+1])

fig, axs = plotref(r, T, rho, p, dlnT, dlnrho, dlnp, s, dsdr,\
    d2lnrho, color='k', xminmax=xminmax)

plt.tight_layout() 
    
axs[0,0].set_title('          ' + dirname_stripped, **csfont)
    
plt.savefig(plotdir + dirname_stripped + '_reference_state.png', dpi=300)
plt.show()
