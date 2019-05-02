import numpy as np
import matplotlib.pyplot as plt
from plot_ref import plotref
from arbitrary_atmosphere import arbitrary_atmosphere

texdir = '/Users/loren/Thesis_Notes/Hydrostatic_ideal_gas_short/'

plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import basic_constants as bc

# Define radius of a hundred grid points
ri = 3.4139791e10
rm = bc.ri
ro = bc.ro
cp = bc.cp

Tm = bc.T_i
pm = bc.p_i
rhom = bc.rho_i
gam = bc.gamma

r = np.linspace(ri, ro, 1000)
g = bc.G*bc.M/r**2
dgdr = -2.0*g/r

# First make a plot for different values of delta

colors = ['b', 'r', 'g', 'm', 'c', 'k', 'y', 'b', 'r', 'g', 'm', 'c', 'k', 'y']

delta = 0.005*rm
kvals = np.linspace(0, 3.6, 10)
count = 0
firstplot = True
for k in kvals:
    d2sdr2 = np.zeros_like(r)
    dsdr = k*cp/rm*(0.5*(1.0 - np.tanh((r - rm)/delta)))
    s = k*cp*(0.5*((r/rm - 1.0) - (delta/rm)*np.log(np.cosh((r - rm)/delta))))
    
    T, rho, p, dlnT, dlnrho, dlnp, d2lnrho =\
        arbitrary_atmosphere(r, s, dsdr, d2sdr2, g, dgdr, rm, Tm, pm, cp, gam)
    
    if firstplot:
        fig, axs = plotref(r, T, rho, p, dlnT, dlnrho, dlnp, s, dsdr,\
            d2lnrho, label=r'$k=%.1f$' %k, color=colors[count])
        firstplot = False
    else:
        plotref(r, T, rho, p, dlnT, dlnrho, dlnp, s, dsdr, d2lnrho, \
                label=r'$k=%.1f$' %k, color=colors[count],\
                fig=fig, axs=axs)
    count += 1

plt.legend()
plt.tight_layout() 
    
axs[0,0].set_title('          Constant sgrad in RZ, del=0.005', **csfont)
    
plt.savefig(texdir + 'const_sgrad_RZ-CZ_k.pdf')
plt.close()