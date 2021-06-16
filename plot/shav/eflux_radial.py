###############################################
# Author: Loren Matilsky
# Date created: 02/14/2018
#
# This script plots the radial energy fluxes as functions of
# radius using from the Shell_Avgs data

import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import *
from compute_grid_info import compute_theta_grid
from rayleigh_diagnostics import GridInfo
from reference_tools import equation_coefficients

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# allowed args + defaults
kwargs_default = {**script_lineplot_kwargs_default}
kwargs_default.update(lineplot_kwargs_default)
kw = update_dict(kwargs_default, clas)
find_bad_keys(kwargs_default, clas, clas0['routinename'], justwarn=True)

# get data
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], 'Shell_Avgs')
print ('Getting data from ' + kw.the_file)
di = get_dict(kw.the_file)
vals = di['vals']
lut = di['lut']
qv = di['qv']
di_grid = get_grid_info(dirname)
rr = di_grid['rr']
nr = di_grid['nr']

eflux = vals[:, 0, lut[1455]]
hflux = vals[:, 0, lut[1433]]
if True in np.isnan(hflux):
    print ("OUTPUT HEAT FLUX (1433, vol_heat_flux) HAS NANs!!!!")
    print ("Computing manually from discrete integral")
    hflux = np.zeros_like(eflux)
    eq = get_eq(dirname)
    rr = eq.radius
    rr2 = rr**2.
    #rho = eq.density
    #temp = eq.temperature
    heat = eq.heating
    for ir in range(eq.nr - 2, -1, -1):
        #mean_rho = 0.5*(rho[ir] + rho[ir+1])
        #mean_t = 0.5*(temp[ir] + temp[ir+1])
        #mean_t = 0.5*(temp[ir] + temp[ir+1])
        mean_r2 = 0.5*(rr2[ir] + rr2[ir+1])
        mean_q = 0.5*(heat[ir] + heat[ir+1])
        mean_dr = rr[ir] - rr[ir+1]
        fpr2dr = 4.*np.pi*mean_r2*mean_dr
        hflux[ir] = hflux[ir+1] + mean_q*fpr2dr
    hflux = (hflux[0] - hflux)/(4.*np.pi*rr2)

cflux = vals[:, 0, lut[1470]]
kflux = vals[:, 0, lut[1923]]
vflux = -vals[:, 0, lut[1935]] # remember the minus sign in vflux
tflux = hflux + eflux + cflux + kflux + vflux # compute the total flux

# break the enthalpy flux into mean and fluctuating
if 1458 in qv:
    print ("getting enthaly flux (pp) from q = 1458")
    eflux_fluc = vals[:, 0, lut[1458]]
    eflux_mean = eflux - eflux_fluc
else: # do the decomposition "by hand"
    eq = get_eq(dirname)
    nt = di_grid['ntheta']
    tt, tw = di_grid['tt'], di_grid['tw'] 
    tw_2d = tw.reshape((nt, 1))
    rho = (eq.density).reshape((1, nr))
    ref_prs = (eq.pressure).reshape((1, nr))
    ref_temp = (eq.temperature).reshape((1, nr))

    print ("getting enthaly flux (pp) indirectly from the AZ_Avgs file")
    print ('Getting data from ' + AZ_Avgs_file)
    AZ_Avgs_file = get_widest_range_file(clas0['datadir'], 'AZ_Avgs')
    di_az = get_dict(AZ_Avgs_file)
    vals_az = di_az['vals']
    lut_az = di_az['lut']

    vr_av = vals_az[:, :, lut_az[1]]
    entropy_av = vals_az[:, :, lut_az[501]]
    prs_av = vals_az[:, :, lut_az[502]]

    # Calculate mean temp. from EOS
    temp_av = ref_temp*((1.-1./thermo_gamma)*(prs_av/ref_prs) +\
            entropy_av/c_P)

    # And, finally, the enthalpy flux from mean/fluc flows
    eflux_mean_az = rho*c_P*vr_av*temp_av
    
    eflux_mean = np.sum(eflux_mean_az*tw_2d, axis=0)
    eflux_fluc = eflux - eflux_mean

if clas0['magnetism']
    # A Space Oddysey is actually (-4*pi) TIMES the correct Poynting flux
    mflux = -vals[:, 0, lut[2001]]/(4*np.pi)
    tflux += mflux

fpr = 4*np.pi*rr**2

# Create the plot; start with plotting all the energy fluxes
lstar = get_lum(dirname)

plt.plot(rr_n, hflux*fpr/lstar, label=r'$\rm{F}_{heat}$', linewidth=lw)
plt.plot(rr_n, eflux*fpr/lstar, 'm', label = r'$\rm{F}_{enth}$',\
        linewidth=lw)
plt.plot(rr_n, cflux*fpr/lstar, label = r'$\rm{F}_{cond}$', linewidth=lw)
plt.plot(rr_n, kflux*fpr/lstar, label = r'$\rm{F}_{KE}$', linewidth=lw)
plt.plot(rr_n, vflux*fpr/lstar, label = r'$\rm{F}_{visc}$', linewidth=lw)
plt.plot(rr_n, tflux*fpr/lstar, label= r'$\rm{F}_{total}$',\
        linewidth=lw, color='black')
if magnetism:
    plt.plot(rr_n, mflux*fpr/lstar, label=r'$\rm{F}_{Poynting}$',\
        linewidth=lw)
plt.plot(rr_n, eflux_fluc*fpr/lstar, 'm--',\
        label=r'$\rm{F}_{enth,\ pp}$', linewidth=lw)
plt.plot(rr_n, eflux_mean*fpr/lstar, 'm:',\
        label=r'$\rm{F}_{enth,\ mm}$', linewidth=lw)
    plt.xlabel(r'$r/R_\odot$',fontsize=kw.fontsize)
else:
    plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=kw.fontsize)

# Try to find the BCZ from where enthalpy flux goes negative, if desired
# avoid the outer boundary
yvals = np.linspace(minmax[0], minmax[1], 100)
if mark_bcz:
    irneg = np.argmin(eflux[20:] > 0) + 20
    irpos = np.argmin(eflux[20:] > 0) + 19
    rrneg = rr_n[irneg]
    rrpos = rr_n[irpos]
    efluxneg = eflux[irneg]
    efluxpos = eflux[irpos]
    slope =  (efluxpos - efluxneg)/(rrpos - rrneg)
    rbcz = rrpos - efluxpos/slope
    plt.plot(rbcz + np.zeros(100), yvals, 'k--', linewidth=lw)

# Mark radii if desired
if not rvals is None:
    for rval in rvals:
        if rnorm is None:
            rval_n = rval/rsun
        else:
            rval_n = rval/rnorm
#        plt.ylim(ymin, ymax)
        plt.plot(rval_n + np.zeros(100), yvals, 'k--', linewidth=lw)


plt.ylabel(r'$4\pi r^2\ \rm{\times \ (energy \ flux)}\ /\ L_*$',\
        fontsize=12, **csfont)

# Label trace interval
if rotation:
    time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.1f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'
else:
    time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.3f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'

# Make title
the_title = dirname_stripped + '\n' + 'radial energy flux, ' + time_string
if mark_bcz:
    the_title += ('\n' + r'$r_{BCZ}/\rm{rnorm} = %.3f$' %rbcz)

plt.title(the_title, **csfont)

# Create a see-through legend
plt.legend(loc='lower left', shadow=True, ncol=3, fontsize=10)

# Last command
plt.tight_layout()

# Save the plot
print ('Saving the eflux plot at ' + plotdir + savename)
plt.savefig(plotdir + savename, dpi=dpi)

# Show the plot
plt.show()
