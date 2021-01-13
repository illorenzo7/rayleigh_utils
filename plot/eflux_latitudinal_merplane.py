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
from azav_util import plot_azav
from common import *

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
plotlatlines = True
plotboundary = True
minmax = None
linthresh = None
linscale = None
minmaxrz = None
linthreshrz = None
linscalerz = None
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
forced = False
rvals = []
rbcz = None
symlog = False

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-minmaxrz':
        minmaxrz = float(args[i+1]), float(args[i+2])
    elif arg == '-rbcz':
        rbcz = float(args[i+1])
    elif arg == '-noshow':
        showplot = False
    elif arg == '-nosave':
        saveplot = False
    elif arg == '-nocontour':
        plotcontours = False
    elif arg == '-usefile':
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
    elif arg == '-forced':
        forced = True
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for rval_str in rvals_str:
            rvals.append(float(rval_str))
    elif arg == '-symlog':
        symlog = True
    elif arg == '-linthresh':
        linthresh = float(args[i+1])
    elif arg == '-linscale':
        linscale = float(args[i+1])
    elif arg == '-linthreshrz':
        linthreshrz = float(args[i+1])
    elif arg == '-linscalerz':
        linscalerz = float(args[i+1])
    elif arg == '-nocontour':
        plotcontours = False
    elif arg == '-nobound':
        plotboundary = False
    elif arg == '-nolat':
        plotlatlines = False

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

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']
nq = di['nq']

# Get the time range in sec
t1 = translate_times(iter1, dirname, translate_from='iter')['val_sec']
t2 = translate_times(iter2, dirname, translate_from='iter')['val_sec']

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'
 
ind_enth = lut[1456] 
ind_cond = lut[1471]
ind_visc = lut[1936] # might get minus sign from error?
ind_ke = lut[1924]

eft_enth = vals[:, :, ind_enth]
eft_cond = vals[:, :, ind_cond]
eft_visc = -vals[:, :, ind_visc]
eft_ke = vals[:, :, ind_ke]

# Check to see if enthalpy flux from the fluctuating flows 
# was already output
if lut[1459] < nq:
    eft_enth_fluc = vals[:, :, lut[1459]]
    eft_enth_mean = eft_enth - eft_enth_fluc
else: # do the Reynolds decomposition "by hand"
    # Compute the enthalpy flux from mean flows (MER. CIRC.)
    eq = get_eq(dirname)
    rho = (eq.density).reshape((1, nr))
    ref_prs = (eq.pressure).reshape((1, nr))
    ref_temp = (eq.temperature).reshape((1, nr))
    prs_spec_heat = 3.5e8
    gamma = 5./3.
    vt_av = vals[:, :, lut[2]]
    entropy_av = vals[:, :, lut[501]]
    prs_av = vals[:, :, lut[502]]

    # Calculate mean temp. from EOS
    temp_av = ref_temp*((1.-1./gamma)*(prs_av/ref_prs) +\
            entropy_av/prs_spec_heat)

    # And, finally, the enthalpy flux from mean/fluc flows
    eft_enth_mean = rho*prs_spec_heat*vt_av*temp_av
    eft_enth_fluc = eft_enth - eft_enth_mean

eft_tot = eft_enth + eft_cond + eft_visc + eft_ke

max_sig = max(np.std(eft_enth), np.std(eft_cond), np.std(eft_visc),\
              np.std(eft_ke))

if magnetism:
    ind_Poynt = lut[2002]
    eft_Poynt = vals[:, :, ind_Poynt]       
    max_sig = max(max_sig, np.std(eft_Poynt))
    eft_tot += eft_Poynt
max_sig = max(max_sig, np.std(eft_tot))
    
if minmax == 'same':
    minmax = -3.*max_sig, 3.*max_sig

# Set up the actual figure from scratch
fig_width_inches = 7. # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1./8. # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_top_inches = 1 # wider top margin to accommodate subplot titles AND metadata
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)
margin_subplot_top_inches = 1/4 # margin to accommodate just subplot titles
nplots = 7 + magnetism
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

eft_terms = [eft_enth_mean, eft_enth_fluc, eft_enth,\
        eft_cond, eft_visc, eft_ke, eft_tot]

titles = [r'$(\mathbf{\mathcal{F}}_{\rm{enth,mm}})_\theta$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{enth,pp}})_\theta$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{enth}})_\theta$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{cond}})_\theta$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{visc}})_\theta$',\
          r'$(\mathbf{\mathcal{F}}_{\rm{KE}})_\theta$',
          r'$(\mathbf{\mathcal{F}}_{\rm{tot}})_\theta$']
units = r'$\rm{g}\ \rm{s}^{-3}$'

if magnetism:
    eft_terms.insert(6, eft_Poynt)
    titles.insert(6, r'$(\mathbf{\mathcal{F}}_{\rm{Poynt}})_\theta$')

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

for iplot in range(nplots):
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1 - margin_top - subplot_height - margin_subplot_top -\
            (iplot//ncol)*(subplot_height + margin_subplot_top +\
            margin_bottom)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    plot_azav (eft_terms[iplot], rr, cost, fig=fig, ax=ax, units=units,\
           minmax=minmax, plotcontours=plotcontours, rvals=rvals,\
           minmaxrz=minmaxrz, rbcz=rbcz, symlog=symlog,\
    linthresh=linthresh, linscale=linscale, linthreshrz=linthreshrz,\
    linscalerz=linscalerz, plotlatlines=plotlatlines, plotboundary=plotboundary)
    ax.set_title(titles[iplot], va='bottom', **csfont)

# Label averaging interval
if rotation:
    time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.1f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'
else:
    time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.3f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'

# Put some metadata in upper left
fsize = 12
fig.text(margin_x, 1 - 0.1*margin_top, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.3*margin_top,\
        'Latitudinal energy flux (zonally averaged)', ha='left', va='top',\
        fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.5*margin_top, time_string,\
         ha='left', va='top', fontsize=fsize, **csfont)

savefile = plotdir + dirname_stripped + '_eflux_latitudinal_merplane_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
if saveplot:
    print ('Saving latitudinal energy fluxes (in the meridional plane) at ' +\
       savefile + ' ...')
    plt.savefig(savefile, dpi=300)
if showplot:
    plt.show()
plt.close()
