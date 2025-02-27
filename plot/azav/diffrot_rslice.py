# Author: Loren Matilsky
# Created: 12/19/2022
#
# Description: Script to plot rotation rate along radial lines 

# Import relevant modules
import numpy as np
import pickle
import matplotlib.pyplot as plt
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

# allowed args + defaults
kw_default = dict({'the_file': None, 'latvals': np.array([0., 15., 30., 45., 60., 75.])})

kw_make_figure_default['margin_top_inches'] += 0.25
kw_default.update(kw_make_figure_default)

kw_lineplot_default['legfrac'] = 1/5
kw_lineplot_default['plotleg'] = True
kw_default.update(kw_lineplot_default)

kw = update_dict(kw_default, clas)
kw_make_figure = update_dict(kw_make_figure_default, clas)
kw_lineplot = update_dict(kw_lineplot_default, clas)

if not kw.xcut is None: # make room for label on right
    kw_make_figure.sub_margin_right_inches = default_margin_xlabel
kw_make_figure.margin_top_inches += default_line_height 
# room for diffrot label
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

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
tt_lat = di_grid['tt_lat']
xx = di_grid['xx']

# rotation rate in inertial frame
eq = get_eq(dirname)
omega = eq.omega0 + vp_av/xx

# DR contrast between 0 and 60 degrees
it0, it60_N, it60_S = np.argmin(np.abs(tt_lat)), np.argmin(np.abs(tt_lat - 60)), np.argmin(np.abs(tt_lat + 60))
delta_omega = omega[it0, 0] - (omega[it60_N, 0] + omega[it60_S, 0])/2

fig, axs, fpar = make_figure(**kw_make_figure)
ax = axs[0, 0]

# Plot rotation vs radius at the desired latitudes
profiles = []
kw_lineplot.labels = []
for latval in kw.latvals:
    ilat_N = np.argmin(np.abs(tt_lat - latval))
    ilat_S = np.argmin(np.abs(tt_lat + latval))
    latitude = (tt_lat[ilat_N] - tt_lat[ilat_S])/2 
    # (this is the actual value we get)
    omega_vsr = (omega[ilat_N, :] + omega[ilat_S, :])/2
    profiles.append(omega_vsr/eq.omega0)
    kw_lineplot.labels.append(r'$\rm{%2.1f}$' %latitude + r'$^\circ$')

# x and y labels
kw_lineplot.xlabel = 'radius'
kw_lineplot.ylabel = 'rotation rate'

# show the frame rotation rate
#kw_lineplot.yvals = make_array(kw_lineplot.yvals, tolist=True)
#kw_lineplot.yvals.append(Om0)
lineplot(rr, profiles, ax, **kw_lineplot)

# make title 
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2) 

maintitle = dirname_stripped + '\n' +  r'$\Omega^*/\Omega_0 - 1$' + ', radial slices\n' + time_string + '\n' + r'$\Delta\Omega_{\rm{60}}/\Omega_0$' + (' = %1.2e' %(delta_omega/eq.omega0))

margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
fig.text(margin_x, 1 - margin_y, maintitle, ha='left', va='top', fontsize=default_titlesize)

# save the figure
plotdir = my_mkdir(clas0['plotdir'])
savefile = plotdir + 'diffrot_rslice' + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
plt.close()
