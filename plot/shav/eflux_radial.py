# Author: Loren Matilsky
# Created: 01/27/2023
# This script plots the radial energy fluxes as functions of radius
# using the Shell_Avgs data

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
# need some other things if ref_type = 5.
reftype = eq.reference_type
if reftype == 5:
    pr_num = 1./eq.constants[5] # 1/c_6
    ra_num = eq.constants[1]*pr_num # c_2 * Pr
    di_num = eq.constants[7]*ra_num # Ra * c_8
    if clas0['magnetism']:
        di_prm = 1./eq.constants[6] # 1/c_7
    factor = ra_num/(di_num*pr_num)

# allowed args + defaults
lineplot_kwargs_default['legfrac'] = 0.3
lineplot_kwargs_default['plotleg'] = True
make_figure_kwargs_default.update(lineplot_fig_dimensions)

kwargs_default = dict({'the_file': None, 'the_file_az': None, 'mark_bcz': False})
kwargs_default.update(make_figure_kwargs_default)

kw = update_dict(kwargs_default, clas)
kw_lineplot = update_dict(lineplot_kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)

if not kw.xcut is None: # make room for label on right
    kw_make_figure.sub_margin_right_inches = default_margin_xlabel
if kw.mark_bcz: # make room for the bcz label
    kw_make_figure.margin_top_inches += default_line_height
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

print (buff_line)
# deal with enthalpy flux first (hard one)
if reftype == 5: # do some special stuff here
    print ("reference_type = 5")
    print ("computing enthalpy flux from separate pressure and entropy fluxes,")
    print("scaled appropriately")

    # total enthalpy flux
    prsflux = -vals[:, 0, lut[1944]]
    entflux = vals[:, 0, lut[1440]]*factor
    eflux = prsflux + entflux

    prsflux_fluc = -vals[:, 0, lut[1947]]
    entflux_fluc = vals[:, 0, lut[1441]]*factor
    eflux_fluc = prsflux_fluc + entflux_fluc
    
    eflux_mean = eflux - eflux_fluc
else:
    eflux = vals[:, 0, lut[1455]]
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

# heat flux also requires some special care, sometimes has Nans
hflux = vals[:, 0, lut[1433]]
if True in np.isnan(hflux):
    print (buff_line)
    print ("OUTPUT HEAT FLUX (1433, vol_heat_flux) HAS NANs!!!!")
    print ("computing heat flux manually from discrete integral")
    hflux = np.zeros(nr)
    rr2 = rr**2.
    for ir in range(nr - 2, -1, -1):
        mean_rr2 = 0.5*(rr2[ir] + rr2[ir+1])
        mean_heat = 0.5*(eq.heat[ir] + eq.heat[ir+1])
        mean_dr = rr[ir] - rr[ir+1]
        fpr2dr = 4.*np.pi*mean_rr2*mean_dr
        hflux[ir] = hflux[ir+1] + mean_heat*fpr2dr
    hflux = (hflux[0] - hflux)/(4.*np.pi*rr2)

# other fluxes are pretty easy
cflux = vals[:, 0, lut[1470]]
kflux = vals[:, 0, lut[1923]]
vflux = -vals[:, 0, lut[1935]] # remember the minus sign in vflux

# some fluxes need multiplication if ref_type = 5
if reftype == 5:
    hflux *= factor
    cflux *= factor

    # deal with possible "flux ratio" parameter
    the_heating_file = None
    for fname in os.listdir(dirname):
        if 'HeatingCooling' in fname:
            the_heating_file = fname
    if not the_heating_file is None:
        # first, shift the heating by the bottom flux
        shift = np.abs(np.min(hflux))*rr[-1]**2/rr**2
        if 'fluxratio' in the_heating_file: # there is a top flux too,
        # shift the heating by that
            shift += np.max(hflux)*rr[0]**2/rr**2
        hflux += shift

tflux = hflux + eflux + cflux + kflux + vflux # compute the total flux


profiles = [hflux, cflux, eflux, eflux_mean, eflux_fluc, kflux, vflux]
kw_lineplot.labels = ['hflux', 'cflux', 'eflux', 'eflux (mm)', 'eflux (pp)', 'kflux', 'vflux']

if clas0['magnetism']:
    # A Space Oddysey is actually (-4*pi) TIMES the correct Poynting flux
    mflux = -vals[:, 0, lut[2001]]/(4*np.pi)
    profiles.append(mflux)
    kw_lineplot.labels.append('mflux')
    tflux += mflux
profiles.append(tflux)
kw_lineplot.labels.append('total')

# integrate and normalize
#lstar = eq.lum
fpr = 4*np.pi*rr**2
lstar = 0.5*((fpr*hflux)[0] + (fpr*hflux)[-1])
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
