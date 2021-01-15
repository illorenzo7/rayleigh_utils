import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from common import *
        get_satvals, saturate_array, get_widest_range_file, get_dict
from plotcommon import axis_range, default_axes_1by1
from rayleigh_diagnostics import Meridional_Slices, ReferenceState
from varprops import texlabels, var_indices, texunits
from azav_util import plot_azav
from get_parameter import get_parameter

# Get command line arguments
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)
radatadir = dirname + '/Meridional_Slices/'
basedir = dirname + '/plots/merslice/'

file_list, int_file_list, nfiles = get_file_lists(radatadir)

minmax = None
posdef = False
symlog = False
logscale = False
iiter = nfiles - 1 # by default plot the last iteration
varname = 'all' # by default plot all the fluid variables
plotlatlines = True

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-ir':
        ir = int(args[i+1])
    elif arg == '-rval':
        rval = float(args[i+1])
    elif arg == '-var':
        varname = args[i+1]
    elif arg == '-symlog':
        symlog = True
    elif arg == '-log':
        logscale = True
    elif arg == '-iter':
        desired_iter = int(args[i+1])
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-sec':
        time = float(args[i+1])
        di_trans = translate_times(time, dirname, translate_from='sec')
        desired_iter = di_trans['val_iter']
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-day':
        time = float(args[i+1])
        di_trans = translate_times(time, dirname, translate_from='day')
        desired_iter = di_trans['val_iter']
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-prot':
        time = float(args[i+1])
        di_trans = translate_times(time, dirname, translate_from='prot')
        desired_iter = di_trans['val_iter']
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-nolon':
        plotlonlines = False

iter_val = int_file_list[iiter]
fname = file_list[iiter]

# Read in desired equatorial slice
print ("Reading " + radatadir + fname + '...')
mer = Meridional_Slices(radatadir + fname, '')
rr = mer.radius
cost = mer.costheta

# Get AZ Avgs data
datadir = dirname + '/data/'
the_file = get_widest_range_file(datadir, 'AZ_Avgs')
di = get_dict(datadir + the_file)
vals_az = di['vals']
lut_az = di['lut']

var_indices_rev = reverse_dict(var_indices)

if varname == 'all': # plot everything
    index_list = mer.qv
else:
    if varname == 't':
        index_list = ['t_index']
    else:
        try:
            index_list = [var_indices[varname]]
        except:
            index_list = [int(varname)]

# Figure dimensions
subplot_width_inches = 2.5
subplot_height_inches = 5.
margin_inches = 1/8
margin_top_inches = 1 # larger top margin to make room for titles

fig_width_inches = subplot_width_inches + 2*margin_inches
fig_height_inches = subplot_height_inches + margin_top_inches + margin_inches

fig_aspect = fig_height_inches/fig_width_inches
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

for var_index in index_list:
    if var_index in [501, 502, 't_index']:
        varname += '_prime'

    try:
        texlabel = texlabels[varname]
    except:
        texlabel = str(var_index).zfill(4)

    plotdir = basedir + varname + '/'
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)
    
    for iphi in range(mer.nphi):
        if var_index != 't_index':
            field = mer.vals[iphi, :, :, mer.lut[var_index], 0]
        else:
            s_az, p_az = vals_az[:, :, lut_az[501]], vals_az[:, :, lut_az[502]]
            field_s = mer.vals[iphi, :, :, mer.lut[501], 0] - s_az
            field_p = mer.vals[iphi, :, :, mer.lut[502], 0] - p_az
            nt, nr = np.shape(field_s)
            ref = ReferenceState(dirname + '/reference', '')
            ref_p = ref.pressure.reshape((1, nr))
            ref_t = ref.temperature.reshape((1, nr))
            cp = get_parameter(dirname, 'pressure_specific_heat')
            gam = 5./3.
            field = ref_t*(field_p/ref_p*(1 - 1/gam) + field_s/cp)

        if var_index in [1, 2, 3]:
            field /= 100.
        elif var_index in [501, 502]:
            field -= vals_az[:, :, lut_az[var_index]]

        # Get default bounds if not specified
        if minmax is None:
            # Cut away the data near the domain boundaries
            minmax = get_satvals(field, posdef=posdef,\
                logscale=logscale, symlog=symlog)

        # Get longitude in degrees
        lon = mer.phi[iphi]*180/np.pi

        savename = 'merslice_' + varname + '_iter' + fname +\
                ('_lon%03.1f' %lon) + '.png'
        print('Plotting merslice: ' + varname + ', iter ' + fname +\
                (', lon %03f' %lon) + ' ...')

        # create axes
        try:
            units = texunits[varname]
        except:
            units = 'cgs'
        fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
        ax = fig.add_axes((margin_x, margin_y, subplot_width, subplot_height))
        plot_azav (field, rr, cost, fig=fig, ax=ax, units=units,\
                minmax=minmax)

        # Make title + label diff. rot. contrast and no. contours
        fsize = 12
        line_spacing_inches = 1/4
        space = line_spacing_inches/fig_height_inches
        fig.text(margin_x, 1 - space, dirname_stripped, ha='left',\
                va='bottom', fontsize=fsize, **csfont)
        fig.text(margin_x, 1 - 2*space,\
                 'iter = ' + fname, ha='left', va='bottom', fontsize=fsize,\
                 **csfont)
        fig.text(margin_x, 1 - 3*space, texlabel +\
                ('      phi = %03.1f' %lon),\
                 ha='left', va='bottom', fontsize=fsize, **csfont)
        savefile = plotdir + savename
        print ('Saving plot at %s ...' %savefile)
        plt.savefig(savefile, dpi=300)
        plt.close()
