# Author: Loren Matilsky
# Created: 01/29/2019
# This script plots the radial energy fluxes in the meridional plane 
# (convective (enthalpy) -- BOTH mean and fluc, conductive,
# radiative heating, kinetic energy,
# viscous, and Poynting flux (if present) 
# Computes the mean enthalpy flux "manually" from <v_r> and <T>
# In the future, just add quantity codes 1458,1459,1460!
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
    if (arg == '-minmax'):
        minmax = float(args[i+1]), float(args[i+2])
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

# See if magnetism is "on"
try:
    magnetism = get_parameter(dirname, 'magnetism')
except:
    magnetism = False # if magnetism wasn't specified, it must be "off"

# Get AZ_Avgs file
print ('Getting radial energy fluxes from ' + datadir + AZ_Avgs_file + ' ...')
di = get_dict(datadir + AZ_Avgs_file)

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']
nq = di['nq']

# Get necessary grid info
rr = di['rr']
cost = di['cost']
sint = di['sint']
tt = di['tt']
tt_lat = di['tt_lat']
nr, nt = len(rr), len(tt)

# Get flux indices
ind_enth = lut[1455] 
ind_cond = lut[1470]
ind_heat = lut[1433]
ind_visc = lut[1935] # might get minus sign from error?
ind_ke = lut[1923]

efr_enth = vals[:, :, ind_enth]
efr_cond = vals[:, :, ind_cond]
efr_heat = vals[:, :, ind_heat]
efr_visc = -vals[:, :, ind_visc]
efr_ke = vals[:, :, ind_ke]

# Check to see if enthalpy flux from the fluctuating flows 
# was already output
if lut[1458] < nq:
    efr_enth_fluc = vals[:, :, lut[1458]]
    efr_enth_mean = efr_enth - efr_enth_fluc
else: # do the Reynolds decomposition "by hand"
    # Compute the enthalpy flux from mean flows (MER. CIRC.)
    ref = ReferenceState(dirname + '/reference', '')
    rho = (ref.density).reshape((1, nr))
    ref_prs = (ref.pressure).reshape((1, nr))
    ref_temp = (ref.temperature).reshape((1, nr))
    prs_spec_heat = get_parameter(dirname, 'pressure_specific_heat')

    gamma = 5./3.
    vr_av = vals[:, :, lut[1]]
    entropy_av = vals[:, :, lut[501]]
    prs_av = vals[:, :, lut[502]]

    # Calculate mean temp. from EOS
    temp_av = ref_temp*((1.-1./gamma)*(prs_av/ref_prs) +\
            entropy_av/prs_spec_heat)

    # And, finally, the enthalpy flux from mean/fluc flows
    efr_enth_mean = rho*prs_spec_heat*vr_av*temp_av
    efr_enth_fluc = efr_enth - efr_enth_mean

efr_tot = efr_enth + efr_cond + efr_heat + efr_visc + efr_ke

max_sig = max(np.std(efr_enth), np.std(efr_cond), np.std(efr_heat),\
        np.std(efr_visc), np.std(efr_ke))

if magnetism:
    ind_Poynt = lut[2001]
    efr_Poynt = vals[:, :, ind_Poynt]       
    max_sig = max(max_sig, np.std(efr_Poynt))
    efr_tot += efr_Poynt
    
if minmax is None: 
    minmax = -3.*max_sig, 3.*max_sig

# Set up the actual figure from scratch
fig_width_inches = 7. # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1./8. # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_top_inches = 2. # wider top margin to accommodate subplot titles AND metadata
margin_subplot_top_inches = 1. # margin to accommodate just subplot titles
nplots = 8 + magnetism
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

efr_terms = [efr_enth_mean, efr_enth_fluc, efr_enth, efr_cond,\
        efr_heat, efr_visc, efr_ke, efr_tot]

titles = [r'$(\mathbf{\mathcal{F}}_{\rm{enth,mm}})_r$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{enth,pp}})_r$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{enth}})_r$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{cond}})_r$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{heat}})_r$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{visc}})_r$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{KE}})_r$',
          r'$(\mathbf{\mathcal{F}}_{\rm{tot}})_r$']
units = r'$\rm{g}\ \rm{s}^{-3}$'

if magnetism:
    # Insert the magnetism plot to just before the last
    # plot (total flux)
    efr_terms.insert(7, efr_Poynt)
    titles.insert(7, r'$(\mathbf{\mathcal{F}}_{\rm{Poynt}})_r$')

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

for iplot in range(nplots):
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1 - margin_top - subplot_height - \
            (iplot//ncol)*(subplot_height + margin_subplot_top)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    plot_azav (efr_terms[iplot], rr, cost, sint, fig=fig, ax=ax,\
           units=units, minmax=minmax,
           plotcontours=plotcontours)

    ax.set_title(titles[iplot], verticalalignment='bottom', **csfont)

# Put some metadata in upper left
fsize = 12
fig.text(margin_x, 1 - 0.1*margin_top, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.3*margin_top, 'Radial energy flux',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.5*margin_top,\
         str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8),\
         ha='left', va='top', fontsize=fsize, **csfont)

savefile = plotdir + dirname_stripped + '_eflux_radial_merplane_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if saveplot:
    print ('Saving radial energy fluxes (in the meridional plane) at ' +\
       savefile + ' ...')
    plt.savefig(savefile, dpi=300)
if showplot:
    plt.show()
plt.close()
