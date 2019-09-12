# Author: Loren Matilsky
# Date created: 03/02/2019
import matplotlib as mpl
from matplotlib import ticker
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['raco'])
from common import get_file_lists, get_widest_range_file, strip_dirname,\
        rsun, get_dict
from plotcommon import axis_range
from get_parameter import get_parameter

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]

# Data and plot directories
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
nosave = False
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)
dirname_stripped = strip_dirname(dirname)

# Find the time/latitude file(s) the data directory. If there are 
# multiple, by default choose the one with widest range in the trace.
# We need both time-latitude regular (for the B field)
# and time-latitude induction (for the induction terms)
time_latitude_induction_file = get_widest_range_file(datadir,\
        'time-latitude_induction')
time_latitude_file = get_widest_range_file(datadir,\
        'time-latitude')

# more defaults
user_specified_minmax = False
user_specified_xminmax = False
tag = ''

desired_rvals = [0.83] # by default, plot time-radius diagram for fields 
    # mid-CZ (units of solar radius)
navg = 11 # by default average over 11 AZ_Avgs files for each time

# Get command-line arguments
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-minmax'):
        user_specified_minmax = True
        my_min_br = float(args[i+1])
        my_max_br = float(args[i+2])
        my_min_induction = float(args[i+3])
        my_max_induction = float(args[i+4])
    elif (arg == '-nosave'):
        nosave = True

    elif (arg == '-usefiles'):
        time_latitude_file = args[i+1]
        time_latitude_file = time_latitude_file.split('/')[-1]
        time_latitude_induction_file = args[i+2]
        time_latitude_induction_file =\
                time_latitude_induction_file.split('/')[-1]
    elif (arg == '-rvals'):
        string_desired_rvals = args[i+1].split()
        if string_desired_rvals == ['all']:
            desired_rvals = 'all'
        else:
            desired_rvals = []
            for j in range(len(string_desired_rvals)):
                desired_rvals.append(float(string_desired_rvals[j]))
    elif (arg == '-navg'):
        navg = int(args[i+1])
        if (navg % 2 == 0):
            print ("Please don't enter even values for navg!")
            print ("Replacing navg = %i with navg = %i" %(navg, navg + 1))
            navg += 1
    elif (arg == '-xminmax'):
        user_specified_xminmax = True
        xmin = float(args[i+1])
        xmax = float(args[i+2])
    elif (arg == '-tag'):
        tag = '_' + args[i+1]

# Read in the time-latitude data (dictionary form)
# In this part, we really must assume that the user saved
# the separate time-latitude (induction) files wisely, so they all
# have the same values for "times" and "lats", depths, etc.
print ('Getting time-latitude trace from ' + datadir +\
       time_latitude_file + ' ...')
di = get_dict(datadir + time_latitude_file)
print ('Getting time-latitude induction trace from ' + datadir +\
       time_latitude_induction_file + ' ...')
di_ind = get_dict(datadir + time_latitude_induction_file)

vals = di['vals']
vals_ind = di_ind['vals']

times = di['times']
iters = di['iters']
rr = di['rr']
ri = di['ri']; ro = di['ro']; shell_depth = ro - ri
tt_lat = di['tt_lat']
rinds = di['rinds'] # radial locations sampled for the trace
ntheta = di['ntheta']
rvals_sampled = rr[rinds]/rsun

qvals = np.array(di['qvals'])
qvals_ind = np.array(di_ind['qvals'])

niter = di['niter']
nr = di['nr']
nrvals = di['ndepths']

iter1 = di['iter1']
iter2 = di['iter2']

# Get global rotation rate, if present
rotation = get_parameter(dirname, 'rotation')
if rotation:
    angular_velocity = get_parameter(dirname, 'angular_velocity')
    Prot = 2*np.pi/angular_velocity
    tnorm = Prot # normalize the time axis by rotation period if applicable
else:
    tnorm = 86400. # normalize by "days"

br_index = np.argmin(np.abs(qvals - 801))
shear_index = np.argmin(np.abs(qvals_ind - 1616))
compr_index = np.argmin(np.abs(qvals_ind - 1617))
adv_index = np.argmin(np.abs(qvals_ind - 1618))
trans_tot_index = np.argmin(np.abs(qvals_ind - 1619))

i_desiredrvals = []
rvals_to_plot = []
if desired_rvals == 'all':
    i_desiredrvals = np.arange(len(rinds))
    rvals_to_plot = rvals_sampled
else:
    i_desiredrvals = []
    rvals_to_plot = []
    for desired_rval in desired_rvals:
        i_desiredrval = np.argmin(np.abs(rvals_sampled - desired_rval))
        i_desiredrvals.append(i_desiredrval)
        rvals_to_plot.append(rvals_sampled[i_desiredrval])

# Get raw traces of br, btheta, bphi
br = vals[:, :, :, br_index]
shear = vals_ind[:, :, :, shear_index]
compr = vals_ind[:, :, :, compr_index]
adv = vals_ind[:, :, :, adv_index]
trans_tot = vals_ind[:, :, :, trans_tot_index]

# Average these traces in time
over2 = navg//2
br_all = np.zeros((niter - navg + 1, ntheta, nrvals))
shear_all = np.zeros((niter - navg + 1, ntheta, nrvals))
compr_all = np.zeros((niter - navg + 1, ntheta, nrvals))
adv_all = np.zeros((niter - navg + 1, ntheta, nrvals))
trans_tot_all = np.zeros((niter - navg + 1, ntheta, nrvals))

for i in range(navg):
    br_all += br[i:niter - navg + 1 + i]
    shear_all += shear[i:niter - navg + 1 + i]
    compr_all += compr[i:niter - navg + 1 + i]
    adv_all += adv[i:niter - navg + 1 + i]
    trans_tot_all += trans_tot[i:niter - navg + 1 + i]
br_all /= navg
trans_tot_all /= navg
shear_all /= navg
compr_all /= navg
adv_all /= navg

times_trace = times[over2:niter - over2]/tnorm

# Make meshgrid of time/radius
# Take into account if user specified xmin, xmax
if user_specified_xminmax:
    it1 = np.argmin(np.abs(times_trace - xmin))
    it2 = np.argmin(np.abs(times_trace - xmax))
    times_trace = times_trace[it1:it2+1]
    br_all = br_all[it1:it2+1]
    trans_tot_all = trans_tot_all[it1:it2+1]
    shear_all = shear_all[it1:it2+1]
    compr_all = compr_all[it1:it2+1]
    adv_all = adv_all[it1:it2+1]
times2, tt_lat2 = np.meshgrid(times_trace, tt_lat, indexing='ij')

# Loop over the desired radii and save plots
for i in range(len(i_desiredrvals)):
    i_desiredrval = i_desiredrvals[i]
    rval_to_plot = rvals_to_plot[i]
    br_trace = br_all[:, :, i_desiredrval]
    trans_tot_trace = trans_tot_all[:, :, i_desiredrval]
    shear_trace = shear_all[:, :, i_desiredrval]
    compr_trace = compr_all[:, :, i_desiredrval]
    adv_trace = adv_all[:, :, i_desiredrval]
  
    # Make appropriate file name to save
    savename = dirname_stripped + '_time-lat_transport_r_mean_' +\
        ('rval%0.3f_' %rval_to_plot) + str(iter1).zfill(8) + '_' +\
        str(iter2).zfill(8) + tag + '.png'

    if not user_specified_minmax:
        std_br = np.std(br_trace)
        my_min_br, my_max_br = -3.*std_br, 3.*std_br

        std_shear = np.std(shear_trace)
        std_trans_tot = np.std(trans_tot_trace)
        std_compr = np.std(compr_trace)
        std_adv = np.std(adv_trace)
        std_max = max(std_shear, std_trans_tot, std_compr, std_adv)
        my_min_induction, my_max_induction = -3.*std_max, 3.*std_max
     
    # Create figure with  5 panels in a row (time-latitude plots of
    # B_r fluc (prime-prime) production: shear, compression, 
    # advection, total induction, and B_r itself
    fig, axs = plt.subplots(5, 1, figsize=(12, 8), sharex=True, sharey=True)
    ax1 = axs[0]; ax2 = axs[1]; ax3 = axs[2]; ax4 = axs[3]; ax5 = axs[4]

    # first plot: evolution of B_r
    im1 = ax1.pcolormesh(times2, tt_lat2, shear_trace,\
            vmin=my_min_induction, vmax=my_max_induction, cmap='RdYlBu_r')
    im2 = ax2.pcolormesh(times2, tt_lat2, compr_trace,\
            vmin=my_min_induction, vmax=my_max_induction, cmap='RdYlBu_r')
    im3 = ax3.pcolormesh(times2, tt_lat2, adv_trace,\
            vmin=my_min_induction, vmax=my_max_induction, cmap='RdYlBu_r')
    im4 = ax4.pcolormesh(times2, tt_lat2, trans_tot_trace,\
            vmin=my_min_induction, vmax=my_max_induction, cmap='RdYlBu_r')
    im5 = ax5.pcolormesh(times2, tt_lat2, br_trace,\
            vmin=my_min_br, vmax=my_max_br, cmap='RdYlBu_r')

    # Put colorbar next to all plots (possibly normalized separately)
    # First make room and then find location of subplots
    plt.subplots_adjust(left=0.1, right=0.85, wspace=0.03, top=0.9)

    # First, shear:
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax1)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x

    cbar_left = ax_xmax + 0.3*(1 - ax_xmax)
    cbar_bottom = ax_ymin
    cbar_width = 0.07*(1 - ax_xmax)
    cbar_height = ax_delta_y
    cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))
    cax.set_title(r'$\rm{G}\ \rm{s}^{-1}$', **csfont)
    plt.colorbar(im1, cax=cax)

    # Next, trans_tot:
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax4)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x

    cbar_left = ax_xmax + 0.3*(1 - ax_xmax)
    cbar_bottom = ax_ymin
    cbar_width = 0.07*(1 - ax_xmax)
    cbar_height = ax_delta_y
    cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))
    cax.set_title(r'$\rm{G}\ \rm{s}^{-1}$', **csfont)
    plt.colorbar(im4, cax=cax)

    # Next, B_r:
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax5)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x

    cbar_left = ax_xmax + 0.3*(1 - ax_xmax)
    cbar_bottom = ax_ymin
    cbar_width = 0.07*(1 - ax_xmax)
    cbar_height = ax_delta_y
    cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))
    cax.set_title(r'$\rm{G}$', **csfont)
    plt.colorbar(im5, cax=cax)

    # Label x (time) axis
    if rotation:
        timeunit = r'$P_{\rm{rot}}$'
    else:
        timeunit = 'days'

    xlabel = 'time (' + timeunit + ')'
    ax5.set_xlabel(xlabel, **csfont)

    # Label y-axis
    ax3.set_ylabel('latitude (deg.)', **csfont)
    ax3.set_yticks(np.arange(-90, 90, 30))

    # Label the plots by induction terms, , , ... B_r
    ax_xmin, ax_xmax = ax1.get_xlim()
    ax_Dx = ax_xmax - ax_xmin
    label1 = r'$[\overline{\mathbf{B}}\cdot\nabla\overline{\mathbf{v}}]_r$'
    ax1.text(ax_xmin + 1.01*ax_Dx, 0., label1, va='center',\
            rotation=270)
    label1 = r'$-\overline{B_r}(\nabla\cdot\overline{\mathbf{v}})$'
    ax2.text(ax_xmin + 1.01*ax_Dx, 0.,  label1, va='center',\
            rotation=270)
    label1 = r'$-[\overline{\mathbf{v}}\cdot\nabla\overline{\mathbf{B}}]_r$'
    ax3.text(ax_xmin + 1.01*ax_Dx, 0., label1, va='center',\
            rotation=270)
    label1 = r'$[\nabla\times(\overline{\mathbf{v}}\times\overline{\mathbf{B}})]_r$'
    ax4.text(ax_xmin + 1.01*ax_Dx, 0.,  label1, va='center',\
            rotation=270)
    ax5.text(ax_xmin + 1.01*ax_Dx, 0.,  r'$\overline{B_r}$', va='center',\
            rotation=270)

    # Put some useful information on the title
    averaging_time = (times[-1] - times[0])/niter*navg/86400.
    title = dirname_stripped + '     ' +\
            (r'$r/R_\odot\ =\ %0.3f$' %rval_to_plot) + '     ' +\
            ('t_avg = %.1f days' %averaging_time) +\
            '     ' + ('time_unit = %.1f days' %(tnorm/86400.))

    ax1.set_title(title, **csfont)
    # Get ticks everywhere
    plt.sca(ax1)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    plt.sca(ax2)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    plt.sca(ax3)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    plt.sca(ax4)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    plt.sca(ax5)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    # Save the plot
    print ('Saving the time-latitude plot at ' + plotdir +\
            savename + ' ...')

    if not nosave:
        plt.savefig(plotdir + savename, dpi=300)

    # Show the plot if only plotting at one latitude
    if len(rvals_to_plot) == 1:
        plt.show()
