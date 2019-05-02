import numpy as np
import matplotlib.pyplot as plt
from plot_ref_and_dref import plotref

texdir = '/Users/loren/Thesis_Notes/Hydrostatic_ideal_gas_short/'

plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import basic_constants as bc

# Define radius of a hundred grid points
rm = bc.ri
Tm = bc.T_i
pm = bc.T_i
rhom = bc.rho_i

ri = 3.4139791e10

r = np.linspace(ri, rm, 100)

# Set up figure axes
fig, axs = plt.subplots(4, 2, figsize= (10,10), sharex=True)
colors = ['b', 'r', 'g', 'm', 'c', 'k', 'y']
    
# Plot a range of polytropic indices
count = 0
for n in np.arange(1.5, 20, 3):
    # Define the constant a
    a = bc.G*bc.M/((n+1.0)*bc.R*Tm*rm)
    zeta = a*rm/r + (1.0 - a)
    
    T = Tm*zeta
    p = pm*zeta**(n + 1.0)
    rho = rhom*zeta**n
    s = bc.cv*(n/bc.n0 - 1.0)*(np.log(r/rm) - np.log(a + (1.0 - a)*r/rm))
    dsdr = (n/bc.n0 - 1.0)*bc.cv/(r + (1.0 - a)*r**2/(a*rm))
    
    dzeta = -a*rm/r**2.0
    dlnzeta = dzeta/zeta
    dlnT = dlnzeta
    dlnrho = n*dlnzeta
    dlnp = (n + 1.0)*dlnzeta
    
    # Plot the polytrope for this n-value
    
    plotref(fig, axs, r, T, rho, p, dlnT, dlnrho, dlnp, s, dsdr,\
            label=r'$n=%.1f$' %n, color=colors[count])
#    plotref(fig, axs, r, T, rho, p, dlnT, dlnrho, dlnp)  
    count += 1

plt.legend()
plt.tight_layout() 
    
axs[0,0].set_title('    Stable polytrope ' + '(radiative zone)', **csfont)
    
plt.savefig(texdir + 'RZ_polytrope.pdf')
plt.close()