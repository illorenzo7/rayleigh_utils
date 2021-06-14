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
kwargs = dict({'the_file': None, 'minmax': None, 'rvals': None, 'latvals': np.array([0., 15., 30., 45., 60., 75.])})
args = sys.argv
clas0, clas = read_clas(args)
kwargs = dotdict(update_kwargs(clas, kwargs))

dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)
the_file = kwargs.the_file
minmax = kwargs.minmax
rvals = kwargs.rvals
latvals = kwargs.latvals

# get data
if the_file is None:
    the_file = get_widest_range_file(clas0['datadir'], 'AZ_Avgs')
print ('Getting data from ' + the_file)
di = get_dict(the_file)
vals = di['vals']
lut = di['lut']
vp_av = vals[:, :, lut[3]]

# Get necessary grid info
di_grid = get_grid_info(dirname)
rr = di_grid['rr']/rsun
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

# figure dimensions
nplots = 1
sub_width_inches = 4.5
sub_height_inches = 3
margin_top_inches = 1 # larger top margin to make room for title
margin_bottom_inches = 1/2
margin_right_inches = 2 # make room for legend

fig, axs, fpar = make_figure(nplots=nplots, sub_width_inches=sub_width_inches, sub_height_inches=sub_height_inches, margin_top_inches=margin_top_inches, margin_left_inches=default_margin_ylabel, margin_right_inches=margin_right_inches, margin_bottom_inches=default_margin_xlabel)
ax = axs[0, 0]

# Plot rotation vs radius at the desired latitudes
count = 0
for latval in latvals:
    ilat_N = np.argmin(np.abs(tt_lat - latval))
    ilat_S = np.argmin(np.abs(tt_lat + latval))
    latitude = (tt_lat[ilat_N] - tt_lat[ilat_S])/2 
    # (this is the actual value we get)
    Om_vs_r = (Om[ilat_N, :] + Om[ilat_S, :])/2
    lineplot(rr, Om_vs_r, ax, label=r'$\rm{%2.1f}$' %latitude + r'$^\circ$', color=color_order[count], minmax=minmax)
    count += 1

# adjust y axis
if minmax is None:
    minmax = ax.get_ylim()
    ax.set_ylim(minmax)

# show the frame rotation rate
lineplot(rr, Om0 + np.zeros_like(rr), ax, color='k', linestyle='--', label=r'$\Omega_0=%.1f$' %Om0, xlabel=r'$r/R_\odot$', ylabel=r'$\Omega/2\pi$' + ' [nHz]', minmax=minmax, xvals=rvals)

# make title 
fontsize = default_titlesize
iter1, iter2 = get_iters_from_file(the_file)
time_string = get_time_string(dirname, iter1, iter2) 
the_title = dirname_stripped + '\n' +  r'$\Omega(r,\theta)$' + '\n' + time_string + '\n' + r'$\Delta\Omega_{\rm{60}}$' + (' = %.1f nHz' %Delta_Om)
margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
fig.text(margin_x, 1 - margin_y, the_title, ha='left', va='top', fontsize=fontsize)
plt.legend(title='latitude', fontsize=fontsize, loc='upper left', bbox_to_anchor=(1.02, 1))

# save the figure
plotdir = my_mkdir(clas0['plotdir'])
savefile = plotdir + 'diffrot_rslice' + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
plt.close()
