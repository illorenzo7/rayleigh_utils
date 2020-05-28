# Author: Loren Matilsky
# Created: 01/29/2019
# This script plots the zonal energy fluxes in the meridional plane 
# (convective (enthalpy), conductive, kinetic energy, viscous,
# and Poynting flux (if present) 
# ...for the Rayleigh run directory indicated by [dirname]. 
# I have a hunch these should all be zero!
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
from azav_util import plot_azav
from common import get_widest_range_file, strip_dirname, get_dict
from get_parameter import get_parameter
from get_eq import get_eq

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
rbcz = None
minmax = None
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-minmax'):
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-rbcz':
        rbcz = float(args[i+1])
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
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
print ('Getting zonal energy fluxes from ' + datadir +\
       AZ_Avgs_file + ' ...')
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

ind_enth = lut[1457] 
# The azimuthal conductive flux had better average to zero ...
ind_visc = lut[1937] # might get minus sign from error?
ind_ke = lut[1925]

efp_enth = vals[:, :, ind_enth]
efp_visc = -vals[:, :, ind_visc]
efp_ke = vals[:, :, ind_ke]

# Check to see if enthalpy flux from the fluctuating flows 
# was already output
if lut[1460] < nq:
    efp_enth_fluc = vals[:, :, lut[1460]]
    efp_enth_mean = efp_enth - efp_enth_fluc
else: # do the Reynolds decomposition "by hand"
    # Compute the enthalpy flux from mean flows (MER. CIRC.)
    eq = get_eq(dirname)
    rho = (eq.density).reshape((1, nr))
    ref_prs = (eq.pressure).reshape((1, nr))
    ref_temp = (eq.temperature).reshape((1, nr))
    prs_spec_heat = 3.5e8

    gamma = 5./3.
    vp_av = vals[:, :, lut[3]]
    entropy_av = vals[:, :, lut[501]]
    prs_av = vals[:, :, lut[502]]

    # Calculate mean temp. from EOS
    temp_av = ref_temp*((1.-1./gamma)*(prs_av/ref_prs) +\
            entropy_av/prs_spec_heat)

    # And, finally, the enthalpy flux from mean/fluc flows
    efp_enth_mean = rho*prs_spec_heat*vp_av*temp_av
    efp_enth_fluc = efp_enth - efp_enth_mean

efp_tot = efp_enth + efp_visc + efp_ke

max_sig = max(np.std(efp_enth), np.std(efp_visc), np.std(efp_ke))

if magnetism:
    ind_Poynt = lut[2003]
    efp_Poynt = vals[:, :, ind_Poynt]       
    max_sig = max(max_sig, np.std(efp_Poynt))
    efp_tot += efp_Poynt
    
if minmax is None:
    minmax = -3.*max_sig, 3.*max_sig

# Set up the actual figure from scratch
fig_width_inches = 7. # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1./8. # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)
margin_top_inches = 1 # wider top margin to accommodate metadata
margin_subplot_top_inches = 1/4 # margin to accommodate just subplot titles
nplots = 6 + magnetism
ncol = 3 # put three plots per row
nrow = np.int(np.ceil(nplots/3))

subplot_width_inches = (fig_width_inches - (ncol + 1)*margin_inches)/ncol
    # Make the subplot width so that ncol subplots fit together side-by-side
    # with margins in between them and at the left and right.
subplot_height_inches = 2*subplot_width_inches # Each subplot should have an
    # aspect ratio of y/x = 2/1 to accommodate meridional planes. 
fig_height_inches = margin_top_inches + nrow*(subplot_height_inches +\
        margin_subplot_top_inches + margin_bottom_inches)
fig_aspect = fig_height_inches/fig_width_inches

# "Margin" in "figure units"; figure units extend from 0 to 1 in BOTH 
# directions, so unitless dimensions of margin will be different in x and y
# to force an equal physical margin
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches
margin_subplot_top = margin_subplot_top_inches/fig_height_inches

# Subplot dimensions in figure units
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

efp_terms = [efp_enth_mean, efp_enth_fluc, efp_enth,\
        efp_visc, efp_ke, efp_tot]

titles = [r'$(\mathbf{\mathcal{F}}_{\rm{enth,mm}})_\phi$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{enth,pp}})_\phi$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{enth}})_\phi$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{visc}})_\phi$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{KE}})_\phi$',
          r'$(\mathbf{\mathcal{F}}_{\rm{tot}})_\phi$']
units = r'$\rm{g}\ \rm{s}^{-3}$'

if magnetism:
    efp_terms.insert(5, efp_Poynt)
    titles.insert(5, r'$(\mathbf{\mathcal{F}}_{\rm{Poynt}})_\phi$')

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

for iplot in range(nplots):
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1 - margin_top - subplot_height - margin_subplot_top -\
            (iplot//ncol)*(subplot_height + margin_subplot_top +\
            margin_bottom)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    plot_azav (efp_terms[iplot], rr, cost, fig=fig, ax=ax,
        units = units, minmax=minmax, plotcontours=plotcontours)

    ax.set_title(titles[iplot], va='bottom', **csfont)

# Put some metadata in upper left
fsize = 12
fig.text(margin_x, 1 - 0.1*margin_top, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.3*margin_top, 'Zonal energy flux (zonally averaged)',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.5*margin_top,\
         str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8),\
         ha='left', va='top', fontsize=fsize, **csfont)

savefile = plotdir + dirname_stripped + '_eflux_zonal_merplane_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if saveplot:
    print ('Saving zonal energy fluxes (in the meridional plane) at ' +\
       savefile + ' ...')
    plt.savefig(savefile, dpi=300)
if showplot:
    plt.show()
plt.close()
