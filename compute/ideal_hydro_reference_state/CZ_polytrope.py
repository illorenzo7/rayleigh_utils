import numpy as np
import matplotlib.pyplot as plt
from plot_ref_and_dref import plotref

texdir = '/Users/loren/Thesis_Notes/Hydrostatic_ideal_gas_short/'

plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import basic_constants as bc

# Define radius of a hundred grid points
r = np.linspace(bc.ri, bc.ro, 100)
# Define the constant a0
a0 = bc.G*bc.M/(bc.cp*bc.T_i*bc.ri)
zeta = a0*bc.ri/r + (1.0 - a0)

T = bc.T_i*zeta
p = bc.p_i*zeta**(bc.n0 + 1.0)
rho = bc.rho_i*zeta**bc.n0

dzeta = -a0*bc.ri/r**2.0
dlnzeta = dzeta/zeta
dlnT = dlnzeta
dlnrho = bc.n0*dlnzeta
dlnp = (bc.n0 + 1.0)*dlnzeta

# Entropy gradients (and entropy) are all zero
s = np.zeros_like(r)
dsdr = np.zeros_like(r)

# Plot this basic polytrope
fig, axs = plt.subplots(4, 2, figsize= (10,10), sharex=True)

plotref(fig, axs, r, T, rho, p, dlnT, dlnrho, dlnp, s, dsdr)

axs[0,0].set_title('    Adiabatic polytrope ' + r'$(N_\rho=3)$', **csfont)

axs[0,1].set_ylim(-8e-10, 0)
axs[1,1].set_ylim(-8e-10, 0)
axs[2,1].set_ylim(-8e-10, 0)

plt.tight_layout()

plt.savefig(texdir + 'CZ_polytrope.pdf')
plt.close()