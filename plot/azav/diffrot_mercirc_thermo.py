# Author: Loren Matilsky
# Created: 12/19/2022
#
# Description: Script to plot rotation-rate contours in meridional plane

import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['raco'] + '/reference_state')
from arbitrary_atmosphere import compute_Di_v
from azav_util import *
from common import *
from plotcommon import *
from cla_util import *
from numbers_util import get_dr_contrast

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname, wrap=True)
rotation = clas0['rotation']

# get desired shells to average over for DR numbers
rvals = clas.rvals
if rvals is None:
    rvals = interpret_rvals(dirname, ['rmin', 'rmax'])
nshells = len(rvals) - 1

# allowed args + defaults
# key unique to this script
kwargs_default = dict({'the_file': None, 'rel': True, 'sub': True, 'nrho': 3, 'beta': 0.759, 'gamma': gamma_ideal, 'verbose': False, 'verbose': False})

# also need make figure kwargs
azav_fig_dimensions['margin_top_inches'] += 0.5*nshells
nlines = get_num_lines(clas0.dirname_label)
azav_fig_dimensions['margin_top_inches'] += (nlines-1)*default_line_height
make_figure_kwargs_default.update(azav_fig_dimensions)

# make sure there is enough space for three plots
kw_make_figure = make_figure_kwargs_default.copy()
kw_make_figure['nrow'] = 1 # change the default
kw_make_figure['ncol'] = 3 # change the default

kwargs_default.update(kw_make_figure)

# and of course need plot_azav kwargs
kwargs_default.update(plot_azav_kwargs_default)

# overwrite defaults, first main kwargs
kw = update_dict(kwargs_default, clas)
kw_plot_azav = update_dict(plot_azav_kwargs_default, clas)
kw_make_figure = update_dict(kw_make_figure, clas)

# check for bad keys
find_bad_keys(kwargs_default, clas, clas0['routinename'], justwarn=True)
if not kw.rcut is None:  
    # need room for two colorbars and line up top stating rcut 
    kw_make_figure.margin_top_inches += 1/4
    kw_make_figure.sub_margin_bottom_inches *= 2

# get data
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], 'AZ_Avgs')
print ('Getting data from ' + kw.the_file)
di = get_dict(kw.the_file)
vals = di['vals']
lut = di['lut']
vp_av = vals[:, :, lut[3]]

# get necessary grid info
ntheta = np.shape(vals)[0]
gi = get_grid_info(dirname, ntheta=ntheta)
rr = gi['rr']
cost = gi['cost']
tt_lat = gi['tt_lat']
xx = gi['xx']

# rotation rate in inertial frame
eq = get_eq(dirname)
omega = eq.omega0 + vp_av/xx

# make plot
fig, axs, fpar = make_figure(**kw_make_figure)
print("ax =", axs)
axs = axs.flatten()
ax = axs[0]

# make diff. rot.

kw_plot_azav.plotlatlines = False 
plot_azav ((omega - eq.omega0)/eq.omega0, rr, cost, fig, ax, **kw_plot_azav)

# make title 
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2, threelines=True) 
maintitle = clas0.dirname_label + '\n' +  r'$\Omega$' + ' (rotation rate)\n' + time_string 

# get DR numbers in spherical shells
# loop over shells and add line of text
for ishell in range(nshells):
    # print shell info first
    r1 = rvals[ishell]
    r2 = rvals[ishell+1]

    # then non-D numbers in shell
    diffrot = get_dr_contrast(dirname, r1, r2, the_file=kw.the_file, verbose=kw.verbose, ntheta=ntheta)
    maintitle += '\n' +\
        ( '(' + flt_fmt + ', ' + flt_fmt + '):') %(r1, r2) +\
        '\n' + (r'$\Delta\Omega\ =\ $' + flt_fmt) %diffrot
        #('(r_1, r_2) = (' + flt_fmt + ', ' + flt_fmt + '):') %(r1, r2)# +\

if not kw.rcut is None:
    maintitle += '\nrcut = %1.3e' %kw.rcut
    
margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
fig.text(margin_x, 1 - margin_y, maintitle,\
         ha='left', va='top', fontsize=default_titlesize)

# make the meridional circulation plot

# get density
rho = eq.rho

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    

# get the meridional circulation flows
vr_av, vt_av = vals[:, :, lut[1]], vals[:, :, lut[2]]

# compute the mass flux
rhovm = rho*np.sqrt(vr_av**2 + vt_av**2)

# compute the streamfunction
psi = streamfunction(rho*vr_av, rho*vt_av, rr, cost)

# make CCW negative and CW positive
rhovm *= np.sign(psi)

# make plot
ax = axs[1]

# plot mass flux
kw_plot_azav.plotcontours = False
plot_azav (rhovm, rr, cost, fig, ax, **kw_plot_azav)

# Plot streamfunction contours
lilbit = 0.01
maxabs = np.max(np.abs(psi))
contourlevels = (-maxabs/2., -maxabs/4., -lilbit*maxabs, 0.,\
        lilbit*maxabs, maxabs/4., maxabs/2.)
kw_plot_azav.plotcontours = kw.plotcontours
kw_plot_azav.plotfield = False
kw_plot_azav.contourlevels = contourlevels
plot_azav (psi, rr, cost, fig, ax, **kw_plot_azav)

# make title 
iter1, iter2 = get_iters_from_file(kw.the_file) 
time_string = get_time_string(dirname, iter1, iter2, threelines=True)
margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
maintitle = dirname_stripped + '\n' + 'mass flux (circulation)\n' + time_string
if not kw.rcut is None:
    maintitle += '\nrcut = %1.3e' %kw.rcut

margin_x = fpar['margin_left'] + fpar['sub_margin_left']
width_skip = fpar['sub_width'] + fpar['sub_margin_left'] + fpar['sub_margin_right']
margin_y = default_margin/fpar['height_inches']
fig.text(margin_x + width_skip, 1 - margin_y, maintitle,\
         ha='left', va='top', fontsize=default_titlesize)

# plot entropy perturbation

# get enetropy deviation (harder than it seems)
entr = vals[:, :, lut[501]]

reftype = eq.reference_type

if reftype == 2:
    c_p = get_parameter(dirname, 'pressure_specific_heat')
    k_s = 1./c_p
elif reftype in [4, 5]: 
    c2 = eq.constants[1] # Ra/Pr or Ro_c^2
    diss = compute_Di_v(kw.gamma, kw.beta, kw.nrho)
    if rotation and reftype == 4: 
        # assume my custom states are nondimensionalized by rotational time
        k = 1.19e-5
    else: # assume nondimensionalized by viscous time
        k = 6.38e-12
    k_s = k_T = c2*k

# calculate relative thermal perturbation
entr_rel = k_s*entr

if kw.sub:
    # compute the spherically averaged thermo. vars
    entr_sph = np.sum(entr*gi.tw_2d, axis=0)
    # subtract the spherical mean from the zonal mean
    entr -= entr_sph.reshape((1, gi.nr))

# plot the entropy
ax = axs[2]
#kw_plot_azav.plotfield = True
kw_plot_azav = update_dict(plot_azav_kwargs_default, clas)
plot_azav (entr, rr, cost, fig, ax, **kw_plot_azav)

# make title 
if kw.rel:
    label = r'$\hat{s}/c_p$'
    if kw.sub:
        titletag = '(relative, sub. sph. mean)'
    else:
        titletag = '(relative, full field)'
else:
    label = r'$\hat{s}$'
    if kw.sub:
        titletag = '(dimensional, sub. sph. mean)'
    else:
        titletag = '(dimensional, full field)'

maintitle = dirname_stripped + '\n' +\
        label + '\n' + titletag +\
        '\n' + time_string

if not kw.rcut is None:
    maintitle += '\nrcut = %1.3e' %kw.rcut

fig.text(margin_x + 2*width_skip, 1 - margin_y, maintitle,\
         ha='left', va='top', fontsize=default_titlesize)

# save the figure
plotdir = my_mkdir(clas0['plotdir'])
pretag = ''
if clas0.prepend:
    pretag = dirname_stripped + '_'

savefile = plotdir + pretag + clas0['routinename'] + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
plt.close()
