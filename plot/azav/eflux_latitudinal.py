# Author: Loren Matilsky
# Created: 01/27/2023
# This script plots the "cone-averaged" colatitudinal 
# energy fluxes as functions of theta.
# Since Rayleigh has no "cone-averages", we use the AZ_Avgs.

import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
from common import *
from plotcommon import *
from cla_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# equation coefficients
eq = get_eq(dirname)

# allowed args + defaults
lineplot_kwargs_default['legfrac'] = 0.3
lineplot_kwargs_default['plotleg'] = True
make_figure_kwargs_default.update(lineplot_kwargs_default)

kwargs_default = dict({'the_file': None})
kwargs_default.update(make_figure_kwargs_default)

kw = update_dict(kwargs_default, clas)
kw_lineplot = update_dict(lineplot_kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)

find_bad_keys(kwargs_default, clas, clas0['routinename'], justwarn=True)

# get data
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], 'Shell_Avgs')
print (buff_line)
print ('running ' + clas0['routinename'])
print (buff_line)
print ('reading ' + kw.the_file)
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
    print (buff_line)
    print ("OUTPUT HEAT FLUX (1433, vol_heat_flux) HAS NANs!!!!")
    print ("computing heat flux manually from discrete integral")
    hflux = np.zeros_like(eflux)
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
print (buff_line)
if 1458 in qv:
    print ("getting enthalpy flux (pp) from qval = 1458")
    eflux_fluc = vals[:, 0, lut[1458]]
    eflux_mean = eflux - eflux_fluc
else: # do the decomposition "by hand"
    print ("getting enthalpy flux (pp) indirectly from AZ_Avgs")

    # more grid info
    nt = di_grid['nt']
    tt, tw = di_grid['tt'], di_grid['tw'] 
    tw_2d = tw.reshape((nt, 1))

    # reference state
    eq = get_eq(dirname)
    rho_ref = (eq.rho).reshape((1, nr))
    prs_ref = (eq.prs).reshape((1, nr))
    tmp_ref = (eq.tmp).reshape((1, nr))

    # read AZ_Avgs data
    if kw.the_file_az is None:
        kw.the_file_az = get_widest_range_file(clas0['datadir'], 'AZ_Avgs')
    print ('reading ' + kw.the_file_az)
    di_az = get_dict(kw.the_file_az)
    vals_az = di_az['vals']
    lut_az = di_az['lut']

    vr_av = vals_az[:, :, lut_az[1]]
    s_av = vals_az[:, :, lut_az[501]]
    prs_av = vals_az[:, :, lut_az[502]]

    # calculate mean tmp from EOS
    tmp_av = tmp_ref*( (1.-1./gamma_ideal)*(prs_av/prs_ref) +\
            s_av/eq.c_p )

    # and, finally, the enthalpy flux from mean/fluc flows
    eflux_mean_az = rho_ref*eq.c_p*vr_av*tmp_av
    eflux_mean = np.sum(eflux_mean_az*tw_2d, axis=0)
    eflux_fluc = eflux - eflux_mean

profiles = [eflux, eflux_mean, eflux_fluc, hflux, cflux, vflux]
kw_lineplot.labels = ['eflux', 'eflux (mm)', 'eflux (pp)', 'hflux', 'cflux', 'vflux']

if clas0['magnetism']:
    # A Space Oddysey is actually (-4*pi) TIMES the correct Poynting flux
    mflux = -vals[:, 0, lut[2001]]/(4*np.pi)
    profiles.append(mflux)
    kw_lineplot.labels.append('mflux')
    tflux += mflux
profiles.append(tflux)
kw_lineplot.labels.append('total')

# integrate and normalize
lstar = eq.lum
fpr = 4*np.pi*rr**2
profiles_int = []
for profile in profiles:
    profiles_int.append(fpr*profile/lstar)

# Create the plot; start with plotting all the energy fluxes
fig, axs, fpar = make_figure(**kw_make_figure)
ax = axs[0,0]

# x and y labels
kw_lineplot.xlabel = 'radius'
kw_lineplot.ylabel = r'$4\pi r^2$' + '(flux)/' + r'$L_*$'

# Try to find the BCZ from where enthalpy flux goes negative, if desired
# avoid the outer boundary
if kw.mark_bcz:
    irneg = np.argmin(eflux[20:] > 0) + 20 # argmin gives first place condition breaks down
    irpos = np.argmin(eflux[20:] > 0) + 19 # remember radial indices are reversed
    rrneg = rr[irneg]
    rrpos = rr[irpos]
    efluxneg = eflux[irneg]
    efluxpos = eflux[irpos]
    slope =  (efluxpos - efluxneg)/(rrpos - rrneg)
    rbcz_est = rrpos - efluxpos/slope
    kw_lineplot.xvals = make_array(kw_lineplot.xvals, tolist=True)
    kw_lineplot.xvals.append(rbcz_est)

    # also mark depth of overshoot
    # see where eflux goes positive
    mineflux = np.min(eflux)
    tol = 0.05*mineflux # tol is negative
    irbcz = np.copy(irneg)
    irpos = np.argmin(eflux[irbcz:] < tol) + irbcz
    irneg = np.argmin(eflux[irbcz:] < tol) + irbcz - 1
    rrneg = rr[irneg]
    rrpos = rr[irpos]
    efluxneg = eflux[irneg]
    efluxpos = eflux[irpos]
    slope =  (efluxpos - efluxneg)/(rrpos - rrneg)
    rov_est = rrneg - efluxneg/slope # remember rrneg is above now
    kw_lineplot.xvals.append(rov_est)

lineplot(rr, profiles_int, ax, **kw_lineplot)

# make title 
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2) 
the_title = dirname_stripped + '\n' +  'radial energy flux' + '\n' + time_string
if kw.mark_bcz:
    the_title += ('\n' + r'$r_{BCZ} = %1.3e$' %rbcz_est)
    the_title += ('\n' + r'$r_{os} = %1.3e$' %rov_est)

margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
fig.text(margin_x, 1 - margin_y, the_title,\
         ha='left', va='top', fontsize=default_titlesize)

# save the figure
plotdir = my_mkdir(clas0['plotdir'])
savefile = plotdir + clas0['routinename'] + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

# print the top and bottom fluxes
print (buff_line)
print ("bottom flux (luminosity in):")
print ("%1.5e" %(fpr*tflux)[-1])
print ("bottom heat flux:")
print ("%1.5e" %(fpr*hflux)[-1])

print (buff_line)
print ("top flux (luminosity out):")
print ("%1.5e" %(fpr*tflux)[0])
print ("top  cond. flux")
print ("%1.5e" %(fpr*cflux)[0])

if clas0['saveplot']:
    print (buff_line)
    print ('saving figure at:')
    print(savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
print (buff_line)
