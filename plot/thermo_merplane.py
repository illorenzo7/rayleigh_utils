# Author: Loren Matilsky
# Created: 01/29/2019
# This script plots the average thermodynamic state in the meridional plane:
# pressure, density, temperature, and entropy
# with the spherically symmetric part subtracted out
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
sys.path.append(os.environ['co'])
sys.path.append(os.environ['pl'])
from azav_util import plot_azav
from rayleigh_diagnostics import ReferenceState
from common import get_widest_range_file, strip_dirname, get_dict
from get_parameter import get_parameter
from binormalized_cbar import MidpointNormalize

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
plotcontours = True
minmax = None
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        my_boundstype = 'manual'
        my_min, my_max = float(args[i+1]), float(args[i+2])
        user_specified_minmax = True
    elif arg == '-noshow':
        showplot = False
    elif arg == '-nosave':
        saveplot = False
    elif arg == '-nocontour':
        plotcontours = False
    elif (arg == '-usefile'):
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]

# Get AZ_Avgs file
Shell_Avgs_file = get_widest_range_file(datadir, 'Shell_Avgs') 

print ('Getting zonally averaged thermo. vars from ' + datadir + AZ_Avgs_file + ' ...')
print ('and the spherically averaged thermo. vars from ' + datadir + Shell_Avgs_file + ' ...')

di = get_dict(datadir + AZ_Avgs_file)
di_sph = np.load(datadir + Shell_Avgs_file)

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
vals_sph = di_sph['vals']
lut = di['lut']
lut_sph = di_sph['lut']

# Get grid info
rr, tt, cost, sint = di['rr'], di['tt'], di['cost'], di['sint'] 
nr, nt = di['nr'], di['nt']

# Compute the thermodynamic variables
ref = ReferenceState(dirname + '/reference', '')
ref_rho = (ref.density).reshape((1, nr))
ref_prs = (ref.pressure).reshape((1, nr))
ref_temp = (ref.temperature).reshape((1, nr))
prs_spec_heat = get_parameter(dirname, 'pressure_specific_heat')

try:
    poly_n = get_parameter(dirname, 'poly_n')
except: # assume by default gamma is 5/3
    poly_n = 1.5

# Compute the zonally averaged thermo. vars
entropy_az = vals[:, :, lut[501]]
prs_az = vals[:, :, lut[502]]

# Calculate mean temp. from EOS
temp_az = ref_temp*(prs_az/ref_prs/(poly_n + 1.) + entropy_az/prs_spec_heat)

# Calculate mean density from Ideal Gas Law
rho_az = ref_rho*(prs_az/ref_prs + temp_az/ref_temp)

# Compute the spherically averaged thermo. vars
entropy_sph = (vals_sph[:, lut_sph[501]]).reshape((1, nr))
prs_sph = (vals_sph[:, lut_sph[502]]).reshape((1, nr))
temp_sph = ref_temp*(prs_sph/ref_prs/(poly_n + 1.) + entropy_sph/prs_spec_heat)
rho_sph = ref_rho*(prs_sph/ref_prs + temp_sph/ref_temp)

# Now subtract the spherical mean from the zonal mean
entropy_fluc = entropy_az - entropy_sph
prs_fluc = prs_az - prs_sph
temp_fluc = temp_az - temp_sph
rho_fluc = rho_az - rho_sph

# Set up the actual figure from scratch
fig_width_inches = 7 # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1/8 # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_top_inches = 2 # wider top margin to accommodate subplot titles AND metadata
margin_subplot_top_inches = 1/2 # margin to accommodate just subplot titles
nplots = 4
ncol = 3 # put three plots per row
nrow = np.int(np.ceil(nplots/3))

subplot_width_inches = (fig_width_inches - (ncol + 1)*margin_inches)/ncol
    # Make the subplot width so that ncol subplots fit together side-by-side
    # with margins in between them and at the left and right.
subplot_height_inches = 2*subplot_width_inches # Each subplot should have an
    # aspect ratio of y/x = 2/1 to accommodate meridional planes. 
fig_height_inches = nrow*subplot_height_inches + margin_top_inches +\
    (nrow - 1)*margin_subplot_top_inches + margin_inches 
    # Room for titles on each row and a regular margin on the bottom
fig_aspect = fig_height_inches/fig_width_inches

# "Margin" in "figure units"; figure units extend from 0 to 1 in BOTH 
# directions, so unitless dimensions of margin will be different in x and y
# to force an equal physical margin
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_subplot_top = margin_subplot_top_inches/fig_height_inches

# Subplot dimensions in figure units
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

thermo_terms = [entropy_fluc, prs_fluc, temp_fluc, rho_fluc]

titles = [r'$\langle S\rangle_{\rm{az}} - \langle S \rangle_{\rm{sph}}$',\
        r'$\langle P\rangle_{\rm{az}} - \langle P \rangle_{\rm{sph}}$',\
        r'$\langle T\rangle_{\rm{az}} - \langle T \rangle_{\rm{sph}}$',\
        r'$\langle \rho\rangle_{\rm{az}} - \langle \rho \rangle_{\rm{sph}}$']
units = [r'$\rm{erg}\ \rm{K}^{-1}\ \rm{g}^{-1}$',\
        r'$\rm{dyn}\ \rm{cm}^{-2}$', r'$\rm{K}$', r'$\rm{g}\ \rm{cm}^{-3}$']

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

for iplot in range(nplots):
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1 - margin_top - subplot_height - \
            (iplot//ncol)*(subplot_height + margin_subplot_top)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    plot_azav (thermo_terms[iplot], rr, cost, sint, fig=fig, ax=ax,\
           units=units[iplot], norm=MidpointNormalize(0.),\
           plotcontours=plotcontours)
    ax.set_title(titles[iplot], va='bottom', **csfont)

# Put some metadata in upper left
fsize = 12
fig.text(margin_x, 1 - 0.1*margin_top, 'Thermodynamic state (zonally averaged)',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.3*margin_top, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.5*margin_top,\
         str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8),\
         ha='left', va='top', fontsize=fsize, **csfont)

savefile = plotdir + dirname_stripped + '_thermo_merplane_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if saveplot:
    print ('Saving thermo. vars (in the meridional plane) at ' +\
            savefile + ' ...')
    plt.savefig(savefile, dpi=300)
if showplot:
    plt.show()
plt.close()
