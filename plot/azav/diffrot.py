# Author: Loren Matilsky
# Created: 12/19/2022
#
# Description: Script to plot rotation-rate contours in meridional plane

import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
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

# get desired shells to average over for DR numbers
rvals = clas.rvals
if rvals is None:
    rvals = interpret_rvals(dirname, ['rmin', 'rmax'])
nshells = len(rvals) - 1

# allowed args + defaults
# key unique to this script
kw_default = dict({'the_file': None, 'verbose': False})

# also need make figure kw
azav_fig_dimensions['margin_top_inches'] += 0.5*nshells
nlines = get_num_lines(clas0.dirname_label)
azav_fig_dimensions['margin_top_inches'] += (nlines-1)*default_line_height
kw_make_figure_default.update(azav_fig_dimensions)

kw_default.update(kw_make_figure_default)

# and of course need plot_azav kw
kw_plot_azav_default['plotlatlines'] = False
kw_default.update(kw_plot_azav_default)

# overwrite defaults, first main kw
kw = update_dict(kw_default, clas)
kw_plot_azav = update_dict(kw_plot_azav_default, clas)
kw_make_figure = update_dict(kw_make_figure_default, clas)

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)
if not kw.rcut is None:  
    # need room for two colorbars and line up top stating rcut 
    kw_make_figure.margin_top_inches += 1/4
    kw_make_figure.sub_margin_bottom_inches *= 2

# need to adjust aspect ratio if halfplane is true
if kw.halfplane:
    kw_make_figure.sub_aspect = 1.

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
di_grid = get_grid_info(dirname, ntheta=ntheta)
rr = di_grid['rr']
cost = di_grid['cost']
tt_lat = di_grid['tt_lat']
xx = di_grid['xx']

# frame rate
eq = get_eq(dirname)

# differential rotation in the rotating frame. 
omega = vp_av/xx

# make plot
fig, axs, fpar = make_figure(**kw_make_figure)
ax = axs[0, 0]

plot_azav (omega/eq.omega0, rr, cost, fig, ax, **kw_plot_azav)

# make title 
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2, threelines=True) 
maintitle = clas0.dirname_label + '\n' +  r'$\Omega^*/\Omega_0$' + ' (rotation rate)\n' + time_string 

# get DR numbers in spherical shells
# loop over shells and add line of text
for ishell in range(nshells):
    # print shell info first
    r1 = rvals[ishell]
    r2 = rvals[ishell+1]
    r1_string = optimal_float_string(r1)
    r2_string = optimal_float_string(r2)

    # then non-D numbers in shell
    diffrot = get_dr_contrast(dirname, r1, r2, the_file=kw.the_file, verbose=kw.verbose, ntheta=ntheta)
    maintitle += '\n' +\
        'from r = ' + r1_string + ' to ' + r2_string + ':' +\
        '\n' + (r'$\Delta\Omega\ =\ $' + flt_fmt) %diffrot
        #('(r_1, r_2) = (' + flt_fmt + ', ' + flt_fmt + '):') %(r1, r2)# +\

if not kw.rcut is None:
    maintitle += '\nrcut = %1.3e' %kw.rcut
    
margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
fig.text(margin_x, 1 - margin_y, maintitle,\
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
