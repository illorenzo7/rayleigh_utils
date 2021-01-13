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
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
from azav_util import plot_azav
from common import *

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

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
    elif arg == '-depths':
        rvals = []
        strings = args[i+1].split()
        for st in strings:
            rval = ro - float(st)*d
            rvals.append(rval)
    elif arg == '-rvals':
        rvals = []
        strings = args[i+1].split()
        for st in strings:
            rval = float(st)*rsun
            rvals.append(rval)
    elif arg == '-rvalscm':
        rvals = []
        strings = args[i+1].split()
        for st in strings:
            rval = float(st)
            rvals.append(rval)
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
    elif arg == '-nolats':
        plotlatlines = False

# Get AZ_Avgs file
Shell_Avgs_file = get_widest_range_file(datadir, 'Shell_Avgs') 

print ('Getting zonally averaged thermo. vars from ' + datadir + AZ_Avgs_file + ' ...')
print ('and the spherically averaged thermo. vars from ' + datadir + Shell_Avgs_file + ' ...')

di = get_dict(datadir + AZ_Avgs_file)
di_sph = get_dict(datadir + Shell_Avgs_file)

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
vals_sph = di_sph['vals']
lut = di['lut']
lut_sph = di_sph['lut']

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

# Get grid info
rr, tt, cost, sint = di['rr'], di['tt'], di['cost'], di['sint'] 
nr, nt = di['nr'], di['nt']

# Compute the thermodynamic variables
prs_spec_heat = get_parameter(dirname, 'pressure_specific_heat')

eq = get_eq(dirname)
ref_rho = (eq.density).reshape((1, nr))
ref_prs = (eq.pressure).reshape((1, nr))
ref_temp = (eq.temperature).reshape((1, nr))

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
rho_az = ref_rho*(prs_az/ref_prs - temp_az/ref_temp)

# Compute the spherically averaged thermo. vars
entropy_sph = (vals_sph[:, lut_sph[501]]).reshape((1, nr))
prs_sph = (vals_sph[:, lut_sph[502]]).reshape((1, nr))
temp_sph = ref_temp*(prs_sph/ref_prs/(poly_n + 1.) + entropy_sph/prs_spec_heat)
rho_sph = ref_rho*(prs_sph/ref_prs - temp_sph/ref_temp)

# Now subtract the spherical mean from the zonal mean
entropy_fluc = entropy_az - entropy_sph
prs_fluc = prs_az - prs_sph
temp_fluc = temp_az - temp_sph
rho_fluc = rho_az - rho_sph

# Now divide out the reference profiles (and cp)
entropy_fluc /= prs_spec_heat
prs_fluc /= ref_prs
temp_fluc /= ref_temp
rho_fluc /= ref_rho


# Set up the actual figure from scratch
fig_width_inches = 7 # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1/8 # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)
margin_top_inches = 1 # wider top margin to accommodate subplot titles AND metadata
margin_subplot_top_inches = 1/4 # margin to accommodate just subplot titles
nplots = 4
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

thermo_terms = [entropy_fluc, prs_fluc, temp_fluc, rho_fluc]

titles = [r'$(\langle S\rangle_{\rm{az}} - \langle S \rangle_{\rm{sph}})/c_{\rm{p}}$',\
        r'$(\langle P\rangle_{\rm{az}} - \langle P \rangle_{\rm{sph}})/\overline{P}$',\
        r'$(\langle T\rangle_{\rm{az}} - \langle T \rangle_{\rm{sph}})/\overline{T}$',\
        r'$(\langle \rho\rangle_{\rm{az}} - \langle \rho \rangle_{\rm{sph}})/\overline{\rho}$']

units = ''

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

fsize = 12
for iplot in range(nplots):
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1 - margin_top - subplot_height - margin_subplot_top -\
            (iplot//ncol)*(subplot_height + margin_subplot_top +\
            margin_bottom)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    plot_azav (thermo_terms[iplot], rr, cost, fig=fig, ax=ax, units=units,\
           minmax=minmax, plotcontours=plotcontours, rvals=rvals,\
           minmaxrz=minmaxrz, rbcz=rbcz, symlog=symlog,\
    linthresh=linthresh, linscale=linscale, linthreshrz=linthreshrz,\
    linscalerz=linscalerz, plotlatlines=plotlatlines)
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
fig.text(margin_x, 1 - 0.1*margin_top, 'Thermodynamic state (zonally averaged)',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.3*margin_top, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.5*margin_top, time_string,\
         ha='left', va='top', fontsize=fsize, **csfont)

savefile = plotdir + dirname_stripped + '_thermo_merplane_nond_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if saveplot:
    print ('Saving thermo. vars (in the meridional plane) at ' +\
            savefile + ' ...')
    plt.savefig(savefile, dpi=300)
if showplot:
    plt.show()
plt.close()
