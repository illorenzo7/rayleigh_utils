# Author: Loren Matilsky
# Created: 01/28/2019
# This script plots the latitudinal energy fluxes in the meridional plane 
# (convective (enthalpy), conductive, kinetic energy, viscous,
# and Poynting flux (if present) 
# ...for the Rayleigh run directory indicated by [dirname]. 
# To use an AZ_Avgs file
# different than the one assosciated with the longest averaging range, use
# -usefile [complete name of desired AZ_Avgs file]
# Saves plot in
# [dirname]_eflux_radial_merplane_[first iter]_[last iter].png

import numpy as np
import pickle
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
from azav_util import plot_quiver
from common import *
from get_parameter import get_parameter
from rayleigh_diagnostics import ReferenceState

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Read command-line arguments (CLAs)
showplot = True
saveplot = True
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
scale = 1
scale_by_mag = True
nsample_t = 20
nsample_r = 10
plot_poles = False

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-scale':
        scale = float(args[i+1])
    elif arg == '-noshow':
        showplot = False
    elif arg == '-nosave':
        saveplot = False
    elif arg == '-usefile':
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
    elif arg == '-unitlength':
        scale_by_mag = False
    elif arg == '-nr':
        nsample_r = int(args[i+1])
    elif arg == '-nt':
        nsample_t = int(args[i+1])
    elif arg == '-poles':
        plot_poles = True

# See if magnetism is "on"
try:
    magnetism = get_parameter(dirname, 'magnetism')
except:
    magnetism = False # if magnetism wasn't specified, it must be "off"

# Get AZ_Avgs file
print ('Getting latitudinal energy fluxes from ' + datadir +\
       AZ_Avgs_file + ' ...')
di = get_dict(datadir + AZ_Avgs_file)

# Get necessary grid info
rr = di['rr']
cost = di['cost']
sint = di['sint']
nr, nt = len(rr), len(cost)
rr_2d = rr.reshape((1, nr))

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']
nq = di['nq']

efr_enth = vals[:, :, lut[1455]]
efr_cond = vals[:, :, lut[1470]]
efr_visc = -vals[:, :, lut[1935]]
efr_ke = vals[:, :, lut[1923]]
efr_heat = vals[:, :, lut[1433]]
efr_tot = efr_enth + efr_cond + efr_visc + efr_ke + efr_heat

# Get luminosity
lum = get_parameter(dirname, 'luminosity')
# Subtract that off from efr_tot
print(np.std(efr_tot))
print(np.std(lum/4/np.pi/rr_2d**2))
efr_tot = efr_tot - lum/4/np.pi/rr_2d**2

eft_enth = vals[:, :, lut[1456]]
eft_cond = vals[:, :, lut[1471]]
eft_visc = -vals[:, :, lut[1936]]
eft_ke = vals[:, :, lut[1924]]
eft_tot = eft_enth + eft_cond + eft_visc + eft_ke

if magnetism:
    efr_Poynt = vals[:, :, lut[2001]]       
    eft_Poynt = vals[:, :, lut[2002]]       
    efr_tot += efr_Poynt
    eft_tot += eft_Poynt

# Create plot
subplot_width_inches = 2.5
subplot_height_inches = 5.
margin_inches = 1/8
margin_top_inches = 3/2 # larger top margin to make room for titles

fig_width_inches = subplot_width_inches + 2*margin_inches
fig_height_inches = subplot_height_inches + margin_top_inches + margin_inches

fig_aspect = fig_height_inches/fig_width_inches
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
ax = fig.add_axes((margin_x, margin_y, subplot_width, subplot_height))

plot_quiver (efr_tot, eft_tot, rr, cost, fig=fig, ax=ax,\
        scale_by_mag=scale_by_mag, nsample_r=nsample_r,\
        nsample_t=nsample_t, plot_poles=plot_poles, scale=scale)

# Make title + label diff. rot. contrast and no. contours
fsize = 12
fig.text(margin_x, 1 - 0.05*margin_top, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.2*margin_top, r'$\Omega - \Omega_0$',\
         ha='left', va='top', fontsize=fsize, **csfont)
savefile = plotdir + dirname_stripped + '_eflux_quiver_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + '.png'
print ('Saving plot at %s ...' %savefile)
plt.savefig(savefile, dpi=300)
plt.show()
