import numpy as np
import matplotlib.pyplot as plt
from plot_ref_and_dref import plotref
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

# First make a plot for different values of delta

# Set up figure axes
fig, axs = plt.subplots(4, 2, figsize= (10,10), sharex=True)
k = 1.0
colors = ['b', 'r', 'g', 'm', 'c', 'k', 'y']

deltas = rm*np.array([0.0005, 0.001, 0.01, 0.03, 0.05, 0.1, 0.15])
count = 0
firstplot = True
for delta in deltas:
    d2sdr2 = np.zeros_like(r)
    dsdr = k*cp/rm*(1.0 - np.tanh((r - rm)/delta))
    s = k*cp*((r/rm - 1.0) - (delta/rm)*np.log(np.cosh((r - rm)/delta)))
    g = bc.G*bc.M/r**2
    dgdr = -2.0*g/r
    
    T, rho, p, dlnT, dlnrho, dlnp, d2lnrho =\
        arbitrary_atmosphere(r, s, dsdr, d2sdr2, g, dgdr, rm, Tm, pm, cp, gam)
    
    if firstplot:
        fig, axs = plotref(r, T, rho, p, dlnT, dlnrho, dlnp, s, dsdr,\
            d2lnrho, label=r'$\delta=%1.1e$' %delta, color=colors[count])
        firstplot = False
    else:
        plotref(r, T, rho, p, dlnT, dlnrho, dlnp, s, dsdr, d2lnrho, \
                label=r'$\delta=%1.1e$' %delta, color=colors[count],\
                fig=fig, axs=axs)
    count += 1

plt.legend()
plt.tight_layout() 
    
axs[0,0].set_title('            Constant sgrad in RZ, k=1.0', **csfont)
    
plt.savefig(texdir + 'const_sgrad_RZ-CZ_delta.pdf')
plt.close()