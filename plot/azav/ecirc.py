# Author: Loren Matilsky
# Created: 05/14/2018
# This script plots the meridional circulation cells the meridional plane 
# for the Rayleigh run directory indicated by [dirname], using the AZ_Avgs
# data. To use an AZ_Avgs file different than the one associated with the 
# longest averaging range, run with option
# -usefile [complete name of desired vavg file]
# Saves plot in
# [dirname]_diffrot_[first iter]_[last iter].npy

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
from azav_util import plot_azav, streamfunction
from common import *

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Get density
eq = get_eq(dirname)
rho = eq.density

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Set defaults
nlevs = 15
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
minmax = None
rbcz = None

# Read in CLAs (if any) to change default variable ranges and other options
plotdir = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-rbcz':
        rbcz = float(args[i+1])
    elif arg == '-nlevs':
        nlevs = int(args[i+1])
    elif arg == '-usefile':
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]

# Read in AZ_Avgs data
print ('Getting data from ' + datadir + AZ_Avgs_file + ' ...')
di = get_dict(datadir + AZ_Avgs_file)

vals = di['vals']
lut = di['lut']
iter1, iter2 = di['iter1'], di['iter2']
rr = di['rr']
tt = di['tt']
cost = di['cost']
sint = di['sint']

eflux_r = vals[:, :, lut[1455]]
cflux_r = vals[:, :, lut[1470]]
kflux_r = vals[:, :, lut[1923]]
vflux_r = -vals[:, :, lut[1935]]
hflux = vals[:, :, lut[1433]]
tflux_r = eflux_r + cflux_r + kflux_r + vflux_r + hflux 
# compute the total radial flux
lum = eq.lum
tflux_r -= lum/(4.*np.pi*rr**2.) 
# Subtract out spherically symmetric bit

eflux_t = vals[:, :, lut[1456]]
cflux_t = vals[:, :, lut[1471]]
kflux_t = vals[:, :, lut[1924]]
vflux_t = -vals[:, :, lut[1936]]
tflux_t = eflux_t + cflux_t + kflux_t + vflux_t # compute the total flux

# Compute the magnitude of the flux
flux_mag = np.sqrt(tflux_r**2. + tflux_t**2.)

# Compute the streamfunction 
psi = streamfunction(tflux_r, tflux_t, rr, cost) 

# Make CCW negative and CW positive 
flux_mag *= np.sign(psi)

# Maximum/minimum mass flux over the meridional plane (excluding polar regions)
it15, it75 = np.argmin(np.abs(tt - 11.*np.pi/12.)),\
    np.argmin(np.abs(tt - np.pi/12.)) # ignore problematic poles 
global_min, global_max = np.min(flux_mag[it15:it75, :]),\
        np.max(flux_mag[it15:it75, :])
maxabs = np.max((np.abs(global_min), np.abs(global_max)))

if minmax is None:
    minmax = -maxabs, maxabs

# Create plot
subplot_width_inches = 2.5
subplot_height_inches = 5.
margin_inches = 1/8
margin_top_inches = 1 # larger top margin to make room for titles
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)

fig_width_inches = subplot_width_inches + 2*margin_inches
fig_height_inches = subplot_height_inches + margin_top_inches +\
        margin_bottom_inches

fig_aspect = fig_height_inches/fig_width_inches
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
ax = fig.add_axes((margin_x, margin_bottom, subplot_width, subplot_height))

# Plot energy flux
plot_azav (flux_mag, rr, cost, fig=fig, ax=ax,\
    units = r'$\rm{erg}\ \rm{cm}^{-2}\ \rm{s}^{-1}$', plotcontours=False,\
    minmax=minmax)

# Plot streamfunction contours
lilbit = 0.01
maxabs = np.max(np.abs(psi))
plot_azav (psi, rr, cost, fig=fig, ax=ax, plotfield=False,\
     nlevs=nlevs)
#plot_azav (psi, rr, cost, sint, fig=fig, ax=ax, plotfield=False,\
#     levels=(-maxabs/2., -maxabs/4.,\
#    -lilbit*maxabs, 0., lilbit*maxabs, maxabs/4., maxabs/2.))

# Make title
fsize = 12
fig.text(margin_x, 1 - 1/8*margin_top, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 3/8*margin_top,\
         r'$|\langle\mathbf{\mathcal{F}}_m\rangle|$',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 5/8*margin_top,\
         str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8),\
         ha='left', va='top', fontsize=fsize, **csfont)

savefile = plotdir + dirname_stripped + '_ecirc_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + '.png'
print ('Saving plot at %s ...' %savefile)
plt.savefig(savefile, dpi=300)
plt.show()
