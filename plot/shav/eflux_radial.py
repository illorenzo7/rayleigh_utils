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
sys.path.append(os.environ['rapl'])
from common import *
from plotcommon import *
from cla_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# allowed args + defaults
kwargs_default = dict({'the_file': None, 'mark_bcz': False})
kwargs_default.update(make_figure_kwargs_default)
kwargs_default.update(lineplot_kwargs_default)
kw = update_dict(kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)
kw_lineplot = update_dict(lineplot_kwargs_default, clas)

print ("kw = ", kw_lineplot)
if not kw.xcut is None: # make room for label on right
    kw_make_figure.sub_margin_right_inches = default_margin_xlabel
if kw.mark_bcz: # make room for the bcz label
    kw_make_figure.margin_top_inches += default_line_height
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
    print (buffer_line)
    print ("OUTPUT HEAT FLUX (1433, vol_heat_flux) HAS NANs!!!!")
    print ("Computing manually from discrete integral")
    print (buffer_line)
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
lstar = get_lum(dirname)
fpr = 4*np.pi*rr**2
profiles_int = []
for profile in profiles:
    profiles_int.append(fpr*profile/lstar)

# Create the plot; start with plotting all the energy fluxes
fig, axs, fpar = make_figure(**kw_make_figure)
ax = axs[0,0]

# x and y labels
kw_lineplot.xlabel = r'$r/R_\odot$'
kw_lineplot.ylabel = r'$4\pi r^2$' + '(flux)/' + r'$L_*$'

# Try to find the BCZ from where enthalpy flux goes negative, if desired
# avoid the outer boundary
if kw.mark_bcz:
    irneg = np.argmin(eflux[20:] > 0) + 20
    irpos = np.argmin(eflux[20:] > 0) + 19
    rrneg = rr[irneg]/rsun
    rrpos = rr[irpos]/rsun
    efluxneg = eflux[irneg]
    efluxpos = eflux[irpos]
    slope =  (efluxpos - efluxneg)/(rrpos - rrneg)
    rbcz_est = rrpos - efluxpos/slope
    kw_lineplot.xvals = make_array(kw_lineplot.xvals, tolist=True)
    kw_lineplot.xvals.append(rbcz_est)

lineplot(rr/rsun, profiles_int, ax, **kw_lineplot)

# make title 
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2) 
the_title = dirname_stripped + '\n' +  'radial energy flux' + '\n' + time_string
if kw.mark_bcz:
    the_title += ('\n' + r'$r_{BCZ}/R_\odot = %.3f$' %rbcz_est)
margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
fig.text(margin_x, 1 - margin_y, the_title,\
         ha='left', va='top', fontsize=default_titlesize)

# save the figure
plotdir = my_mkdir(clas0['plotdir'])
savefile = plotdir + clas0['routinename'] + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
