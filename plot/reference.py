#!/Users/loren/anaconda3/bin/python
from get_parameter import *
import numpy as np
import matplotlib.pyplot as plt
import os, sys
from diagnostic_reading import ReferenceState

# Makes use of the ReferenceState class, which has the following attributes:

# self.n_r         : number of radial points
# self.radius      : radial coordinates
# self.density     : density
# self.dlnrho      : logarithmic derivative of density
# self.d2lnrho     : d_by_dr of dlnrho
# self.pressure    : pressure
# self.temperature : temperature
# self.dlnt        : logarithmic derivative of temperature
# self.dsdr        : entropy gradient (radial)
# self.entropy     : entropy
# self.gravity     : gravity

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]

# Directory with data
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'



# Read in reference state and corresponding variables
ref = ReferenceState(dirname + '/reference', '')

rr = ref.radius
rsun = 6.96e10
rr_n = rr/rsun # rr "normalized" by the solar radius
p_ref = ref.pressure
s_ref = ref.entropy
t_ref = ref.temperature
rho_ref = ref.density
Hrho_ref = -1/(ref.dlnrho)/1.e8 # conver cm --> Mm
g_ref = ref.gravity/100. # convert to m/s^2

ro = np.max(rr)
ri = np.min(rr)
d = ro - ri
ir_c = np.argmin(np.abs(rr - (ri + d/2))) #index of the center of the layer

T_c = t_ref[ir_c]
P_c = p_ref[ir_c]
rho_c = rho_ref[ir_c]
S_c = s_ref[ir_c]

plt.plot(rr_n, rho_ref/rho_c, label = (r'$\rho/\rho_c, \rho_c = %.1e$' %rho_c) + ' g/cm' + r'$^3$')
plt.plot(rr_n, t_ref/T_c, label = (r'$T/T_c, T_c = %.1e$' %T_c) + ' K')
plt.plot(rr_n, p_ref/P_c, label = (r'$P/P_c, P_c = %.1e$' %P_c) + ' dyn/cm' + r'$^2$')
plt.plot(rr_n, s_ref/S_c, label = (r'$S/S_c, S_c = %.1e$' %S_c) + ' erg/K/g')

plt.xlabel(r'$r/R_\odot$', fontsize=14)
plt.xlim(ri/rsun, ro/rsun)
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
plt.legend()
plt.tight_layout()
plt.savefig(plotdir + 'reference_td.png', dpi=300)
plt.close()

# Now plot each relevant quantity individually, in its own units (including
# scale height and gravity) 

plt.plot(rr_n, rho_ref)
plt.xlabel(r'$r/R_\odot$', fontsize=14)
plt.ylabel(r'$\hat{\rho}$' + ' (g/cm'  + r'$^3$' + ')', fontsize=14)
plt.xlim(ri/rsun, ro/rsun)
plt.ylim(0, 1.1*np.max(rho_ref))
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
plt.tight_layout()
plt.savefig(plotdir + 'reference_rho.png', dpi=300)
plt.close()

plt.plot(rr_n, p_ref)
plt.xlabel(r'$r/R_\odot$', fontsize=14)
plt.ylabel(r'$\hat{P}$' + ' (dyn/cm'  + r'$^2$' + ')', fontsize=14)
plt.xlim(ri/rsun, ro/rsun)
plt.ylim(0, 1.1*np.max(p_ref))
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
plt.tight_layout()
plt.savefig(plotdir + 'reference_p.png', dpi=300)
plt.close()

plt.plot(rr_n, t_ref)
plt.xlabel(r'$r/R_\odot$', fontsize=14)
plt.ylabel(r'$\hat{T}$' + ' (K)', fontsize=14)
plt.xlim(ri/rsun, ro/rsun)
plt.ylim(0, 1.1*np.max(t_ref))
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
plt.ticklabel_format(scilimits=(0,0))
plt.tight_layout()
plt.savefig(plotdir + 'reference_t.png', dpi=300)
plt.close()

plt.plot(rr_n, s_ref)
plt.xlabel(r'$r/R_\odot$', fontsize=14)
plt.ylabel(r'$\hat{S}$' + ' (erg/K/g)', fontsize=14)
plt.xlim(ri/rsun, ro/rsun)
plt.ylim(0, 1.1*np.max(s_ref))
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
plt.tight_layout()
plt.savefig(plotdir + 'reference_s.png', dpi=300)
plt.close()

plt.plot(rr_n, g_ref)
plt.xlabel(r'$r/R_\odot$', fontsize=14)
plt.ylabel(r'$g(r)$' + ' (m/s' + r'$^2$' + ')', fontsize=14)
plt.xlim(ri/rsun, ro/rsun)
plt.ylim(0, 1.1*np.max(g_ref))
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
plt.tight_layout()
plt.savefig(plotdir + 'reference_g.png', dpi=300)
plt.close()

plt.plot(rr_n, Hrho_ref)
plt.xlabel(r'$r/R_\odot$', fontsize=14)
plt.ylabel(r'$H_\rho(r)$' + ' (Mm)', fontsize=14)
plt.xlim(ri/rsun, ro/rsun)
plt.ylim(0, 1.1*np.max(Hrho_ref))
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
plt.tight_layout()
plt.savefig(plotdir + 'reference_Hrho.png', dpi=300)
plt.close()
