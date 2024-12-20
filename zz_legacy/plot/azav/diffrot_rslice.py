# Author: Loren Matilsky
# Created: 05/14/2018
# This script generates differential rotation plotted along radial lines for
# the Rayleigh run directory indicated by [dirname]. To use  time-averaged 
# AZ_Avgs file different than the one associated with the longest averaging 
# range, use -usefile [complete name of desired vavg file]
# Saves plot in
# [dirname]_diffrot_rslice_[first iter]_[last iter].npy

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
kwargs_default = dict({'the_file': None, 'latvals': np.array([0., 15., 30., 45., 60., 75.])})
kwargs_default.update(make_figure_kwargs_default)
lineplot_kwargs_default['legfrac'] = 1/5
lineplot_kwargs_default['plotleg'] = True
kwargs_default.update(lineplot_kwargs_default)
kw = update_dict(kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)
kw_lineplot = update_dict(lineplot_kwargs_default, clas)

if not kw.xcut is None: # make room for label on right
    kw_make_figure.sub_margin_right_inches = default_margin_xlabel
kw_make_figure.margin_top_inches += default_line_height 
# room for diffrot label
find_bad_keys(kwargs_default, clas, clas0['routinename'], justwarn=True)

# get data
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], 'AZ_Avgs')
print ('Getting data from ' + kw.the_file)
di = get_dict(kw.the_file)
vals = di['vals']
lut = di['lut']
vp_av = vals[:, :, lut[3]]

# Get necessary grid info
di_grid = get_grid_info(dirname)
rr = di_grid['rr']
tt_lat = di_grid['tt_lat']
xx = di_grid['xx']

# Get frame rate rotation and compute differential rotation in the 
# lab frame. 
Om0 = 2*np.pi/compute_Prot(dirname)
Om = vp_av/xx + Om0
Om *= 1e9/2/np.pi # convert from rad/s --> nHz
Om0 *= 1e9/2/np.pi # convert from rad/s --> nHz

# DR contrast between 0 and 60 degrees
it0, it60_N, it60_S = np.argmin(np.abs(tt_lat)), np.argmin(np.abs(tt_lat - 60)), np.argmin(np.abs(tt_lat + 60))
Delta_Om = Om[it0, 0] - (Om[it60_N, 0] + Om[it60_S, 0])/2

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
    Om_vs_r = (Om[ilat_N, :] + Om[ilat_S, :])/2
    profiles.append(Om_vs_r)
    kw_lineplot.labels.append(r'$\rm{%2.1f}$' %latitude + r'$^\circ$')
   
# show the frame rotation rate
kw_lineplot.yvals = make_array(kw_lineplot.yvals, tolist=True)
kw_lineplot.yvals.append(Om0)
lineplot(rr/rsun, profiles, ax, **kw_lineplot)

# make title 
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2) 
the_title = dirname_stripped + '\n' +  r'$\Omega(r,\theta)$' + '\n' + time_string + '\n' + r'$\Delta\Omega_{\rm{60}}$' + (' = %.1f nHz' %Delta_Om)
margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
fig.text(margin_x, 1 - margin_y, the_title, ha='left', va='top', fontsize=default_titlesize)

# save the figure
plotdir = my_mkdir(clas0['plotdir'])
savefile = plotdir + 'diffrot_rslice' + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
plt.close()
