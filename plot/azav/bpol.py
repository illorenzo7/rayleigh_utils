# Author: Loren Matilsky
# Created: 12/19/2022
#
# Description: Script to plot mass flux (stream lines and magnitude)
 
import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from azav_util import *
from common import *
from plotcommon import *
from cla_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname, wrap=True)

# allowed args + defaults
# key unique to this script
kwargs_default = dict({'the_file': None})

# also need make figure kwargs
make_figure_kwargs_default.update(azav_fig_dimensions)
kwargs_default.update(make_figure_kwargs_default)

# of course, plot_azav kwargs, but need to change a few
kwargs_default.update(plot_azav_kwargs_default)

# overwrite defaults
kw = update_dict(kwargs_default, clas)
kw_plot_azav = update_dict(plot_azav_kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)
kw_make_figure.ncol = 3 # room for <B_r> and <B_theta>
kw_make_figure.margin_top_inches += 1/4 # room for letter labels

# check for bad keys
find_bad_keys(kwargs_default, clas, clas0['routinename'], justwarn=True)
if not kw.rcut is None:  
    # need room for two colorbars and line up top stating rcut 
    kw_make_figure.margin_top_inches += 1/4
    kw_make_figure.sub_margin_bottom_inches *= 2

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    

# get data
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], 'AZ_Avgs')
print ('Getting data from ' + kw.the_file)
di = get_dict(kw.the_file)
vals = di['vals']
lut = di['lut']

# Get necessary grid info
ntheta = np.shape(vals)[0]
di_grid = get_grid_info(dirname, ntheta=ntheta)
rr = di_grid['rr']
cost = di_grid['cost']

# get the meridional mean fields
br_av, bt_av = vals[:, :, lut[801]], vals[:, :, lut[802]]

# compute field strength
bpol = np.sqrt(br_av**2 + bt_av**2)

# Compute the streamfunction
psi = streamfunction(br_av, bt_av, rr, cost)

# Make CCW negative and CW positive
bpol *= np.sign(psi)

# make plot
# field strength
fig, axs, fpar = make_figure(**kw_make_figure)
ax = axs[0, 0]
ax.set_title('(a) field lines + strength', fontsize=kw.fontsize, loc='left')
kw_orig = dotdict(kw_plot_azav)
kw_plot_azav.plotcontours = False
plot_azav (bpol, rr, cost, fig, ax, **kw_plot_azav)

# plot field-lines
#lilbit = 0.01
#maxabs = np.max(np.abs(psi))
#contourlevels = (-maxabs/2., -maxabs/4., -lilbit*maxabs, 0.,\
kw_plot_azav.plotcontours = kw.plotcontours
kw_plot_azav.plotfield = False
#kw_plot_azav.contourlevels = contourlevels
plot_azav (psi, rr, cost, fig, ax, **kw_plot_azav)

# plot <B_r> and <B_theta>
kw_plot_azav = kw_orig
kw_plot_azav.plotcontours = False
ax = axs[0, 1]
ax.set_title('(b) ' + r'$\langle B_r\rangle$',\
        fontsize=kw.fontsize, loc='left')
plot_azav (br_av, rr, cost, fig, ax, **kw_plot_azav)

ax = axs[0, 2]
ax.set_title('(c) ' + r'$\langle B_\theta\rangle$',\
        fontsize=kw.fontsize, loc='left')
plot_azav (bt_av, rr, cost, fig, ax, **kw_plot_azav)

# make title 
iter1, iter2 = get_iters_from_file(kw.the_file) 
time_string = get_time_string(dirname, iter1, iter2, threelines=True)
margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
maintitle = dirname_stripped + '\n' + 'poloidal magnetic field\n' + time_string
if not kw.rcut is None:
    maintitle += '\nrcut = %1.3e' %kw.rcut

margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
fig.text(margin_x, 1 - margin_y, maintitle,\
         ha='left', va='top', fontsize=default_titlesize)

# save the figure
plotdir = my_mkdir(clas0['plotdir'])
savefile = plotdir + clas0['routinename'] + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
plt.close()
