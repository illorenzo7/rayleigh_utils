# Author: Loren Matilsky
# Date created: 04/05/2020
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
from time_scales import compute_Prot, compute_tdt
from get_eq import get_eq
from tl_util import plot_tl

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]

# Data and plot directories
datadir = dirname + '/data/'
plotdir = dirname + '/plots/time-lat/'
nosave = False
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)
dirname_stripped = strip_dirname(dirname)

# Find the time/latitude file(s) the data directory. If there are 
# multiple, by default choose the one with widest range in the trace.
the_file = get_widest_range_file(datadir, 'time-latitude_torques')

# more defaults
minmax = None
minmax_Lz = None
xminmax = None
tag = ''

desired_rvals = [0.692] # by default, plot time-radius diagram for fields 
    # mid-RZ (units of solar radius)
navg = 11 # by default average over 11 AZ_Avgs files for each time
forced = False

# Get command-line arguments
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-nosave':
        nosave = True
    elif arg == '-usefile':
        the_file = args[i+1]
        the_file = the_file.split('/')[-1]
    elif arg == '-rvals':
        string_desired_rvals = args[i+1].split()
        if string_desired_rvals == ['all']:
            desired_rvals = 'all'
        else:
            desired_rvals = []
            for j in range(len(string_desired_rvals)):
                desired_rvals.append(float(string_desired_rvals[j]))
    elif arg == '-navg':
        navg = int(args[i+1])
        if (navg % 2 == 0):
            print ("Please don't enter even values for navg!")
            print ("Replacing navg = %i with navg = %i" %(navg, navg + 1))
            navg += 1
    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-lzminmax':
        minmax_Lz = float(args[i+1]), float(args[i+2])
    elif arg == '-tag':
        tag = '_' + args[i+1]
    elif arg == '-forced':
        forced = True

# Read in the time-latitude data (dictionary form)
print ('Getting time-latitude torques trace from ' + datadir +\
       the_file + ' ...')
di = get_dict(datadir + the_file)

vals = di['vals']

times = di['times']
iters = di['iters']
rr = di['rr']
ri = di['ri']; ro = di['ro']; shell_depth = ro - ri
tt_lat = di['tt_lat']
rinds = di['rinds'] # radial locations sampled for the trace
ntheta = di['ntheta']
rvals_sampled = rr[rinds]/rsun

qvals = np.array(di['qvals'])

niter = di['niter']
nr = di['nr']
nrvals = di['ndepths']
nt = di['ntheta']

iter1 = di['iter1']
iter2 = di['iter2']

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

torque_rs_index = np.argmin(np.abs(qvals - 1801))
torque_mm_index = np.argmin(np.abs(qvals - 1802))
torque_cor_index = np.argmin(np.abs(qvals - 1803))
torque_v_index = np.argmin(np.abs(qvals - 1804))
Lz_index = np.argmin(np.abs(qvals - 1819))

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

# Get raw traces of torques and Lz 
torque_mc = vals[:, :, :, torque_mm_index] + vals[:, :, :, torque_cor_index]
torque_rs = vals[:, :, :, torque_rs_index] 
torque_v = vals[:, :, :, torque_v_index] 
Lz = vals[:, :, :, Lz_index] 
if forced:
    vp_index = np.argmin(np.abs(qvals - 3))
    vp = vals[:, :, :, vp_index]
    tacho_r = get_parameter(dirname, 'tacho_r')
    tacho_dr = get_parameter(dirname, 'tacho_dr')
    tacho_tau = get_parameter(dirname, 'tacho_tau')
    forcing = np.zeros((niter, nt, nrvals))
    eq = get_eq(dirname)
    rho = eq.rho

    if os.path.exists(dirname + '/inner_vp'):
        print ("inner_vp file exists, so I assume you have a forcing function which\n quartically matches on to a CZ differential rotation")
        inner_vp = read_inner_vp(dirname + '/inner_vp')
        for it in range(nt):
            for ir in rinds:
                if rr[ir] <= tacho_r:
                    desired_vp = inner_vp[it]*(1.0 - ( (rr[ir] - tacho_r)/(tacho_dr*rr[0]) )**2)**2
                else:
                    desired_vp = 0.0
                forcing[:, it, ir] = -rho[ir]*(vp[:, it, ir] - desired_vp)/tacho_tau
    else:
        forcing_coeff = -rho/tacho_tau*0.5*(1.0 - np.tanh((rr - tacho_r)/(tacho_dr*rr[0])))
        forcing = forcing_coeff[rinds].reshape((1, 1, nrvals))*vp

    # convert forcing function into a torque
    torque_forcing = di['xx'][:, rinds].reshape((1, nt, nrvals))*forcing

# Average these traces in time
# "_all" means "all depths"
over2 = navg//2
torque_rs_all = np.zeros((niter - navg + 1, ntheta, nrvals))
torque_mc_all = np.zeros((niter - navg + 1, ntheta, nrvals))
torque_v_all = np.zeros((niter - navg + 1, ntheta, nrvals))
Lz_all = np.zeros((niter - navg + 1, ntheta, nrvals))
if forced:
    torque_forcing_all = np.zeros((niter - navg + 1, ntheta, nrvals))

for i in range(navg):
    torque_rs_all += torque_rs[i:niter - navg + 1 + i]
    torque_mc_all += torque_mc[i:niter - navg + 1 + i]
    torque_v_all += torque_v[i:niter - navg + 1 + i]
    Lz_all += Lz[i:niter - navg + 1 + i]
    if forced:
        torque_forcing_all += torque_forcing[i:niter - navg + 1 + i]
torque_rs_all /= navg
torque_mc_all /= navg
torque_v_all /= navg
Lz_all /= navg
if forced:
    torque_forcing_all /= navg

times = times[over2:niter - over2]/time_unit

# Make meshgrid of time/radius
# Take into account if user specified xmin, xmax
if xminmax is None: # By default use all times available
    it1 = 0
    it2 = niter - navg
else:
    it1 = np.argmin(np.abs(times - xminmax[0]))
    it2 = np.argmin(np.abs(times - xminmax[1]))
times = times[it1:it2+1]
torque_rs_all = torque_rs_all[it1:it2+1]
torque_mc_all = torque_mc_all[it1:it2+1]
torque_v_all = torque_v_all[it1:it2+1]
Lz_all = Lz_all[it1:it2+1]
if forced:
    torque_forcing_all = torque_forcing_all[it1:it2+1]

times_2d, tt_lat_2d = np.meshgrid(times, tt_lat, indexing='ij')

# Loop over the desired radii and save plots
for i in range(len(i_desiredrvals)):
    i_desiredrval = i_desiredrvals[i]
    rval_to_plot = rvals_to_plot[i]
    torque_rs_loc = torque_rs_all[:, :, i_desiredrval]
    torque_mc_loc = torque_mc_all[:, :, i_desiredrval]
    torque_v_loc = torque_v_all[:, :, i_desiredrval]
    torque_tot_loc = torque_rs_loc + torque_mc_loc + torque_v_loc
    Lz_loc = Lz_all[:, :, i_desiredrval]

    if forced:
        torque_forcing_loc = torque_forcing_all[:, :, i_desiredrval]
        torque_tot_loc += torque_forcing_loc
 
    # Make appropriate file name to save
    savename = dirname_stripped + '_time-lat_torques_' +\
        ('rval%0.3f_' %rval_to_plot) + str(iter1).zfill(8) + '_' +\
        str(iter2).zfill(8) + tag + '.png'

    # Create figure with  5-6 panels in a row 
    # (torques: rs, mc, v, maybe forcing, tot, and L_z)
    nrow = 5 + forced

    # set figure dimensions
    fig_width_inches = 7. + 1./4.
    margin_inches = 1./4.
    margin_bottom_inches = 1./2. # space for x-axis label
    margin_top_inches = 1./2.
    margin_left_inches = 5./8. # space for latitude label
    margin_right_inches = 1.
    subplot_width_inches = fig_width_inches - margin_right_inches -\
            margin_left_inches
    subplot_height_inches = 1.5

    fig_height_inches = nrow*subplot_height_inches +\
            (nrow - 1)*margin_inches + margin_bottom_inches +\
            margin_top_inches

    margin_x = margin_inches/fig_width_inches
    margin_y = margin_inches/fig_height_inches
    subplot_width = subplot_width_inches/fig_width_inches
    subplot_height = subplot_height_inches/fig_height_inches
    margin_left = margin_left_inches/fig_width_inches
    margin_bottom = margin_bottom_inches/fig_height_inches
    margin_top = margin_top_inches/fig_height_inches

    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    
    fields = [torque_rs_loc, torque_mc_loc, torque_v_loc, torque_tot_loc,\
            Lz_loc]
    units = [r'$\rm{g\ cm^{-1}\ s^{-2}}$',\
        r'$\rm{g\ cm^{-1}\ s^{-2}}$',\
        r'$\rm{g\ cm^{-1}\ s^{-2}}$',\
        r'$\rm{g\ cm^{-1}\ s^{-2}}$',\
        r'$\rm{g\ cm^{-1}\ s^{-1}}$']

    if forced:
        fields.insert(2, torque_forcing_loc)
        units.insert(2, r'$\rm{g\ cm^{-1}\ s^{-2}}$')
        mins_and_maxes.insert(2, minmax)

    axs = []
    for irow in range(nrow): 
        if irow < nrow - 1: 
            minmax_loc = minmax
        else:
            minmax_loc = minmax_Lz
        axs.append(fig.add_axes((margin_left, 1. - margin_top -\
                subplot_height - irow*(subplot_height + margin_y),\
                subplot_width, subplot_height)))
        plot_tl(fields[irow], times, tt_lat, minmax=minmax_loc,\
                xminmax=xminmax, navg=navg, fig=fig, ax=axs[irow],\
                units=units[irow])
        if irow < nrow - 1:
            axs[irow].set_xticklabels([])

    # Put colorbar next to all plots (possibly normalized separately)
    # First make room and then find location of subplots
#    plt.subplots_adjust(left=0.1, right=0.85, wspace=0.03, top=0.9)


    # Label x (time) axis
    xlabel = 'time (' + time_label + ')'
    axs[-1].set_xlabel(xlabel, **csfont)

    # Label y-axis
    axs[2].set_ylabel('latitude (deg.)', **csfont)
    axs[2].set_yticks(np.arange(-90, 90, 30))

    # Label the plots by torques and L_z
    labels = [r'$\tau_{\rm{rs}}$', r'$\tau_{\rm{mc}}$',\
            r'$\tau_{\rm{v}}$', r'$\tau_{\rm{tot}}$', r'$\mathcal{L}_z$']
    if forced:
        labels.insert(2, r'$\tau_{\rm{forcing}}$')

    for irow in range(nrow):
        label = labels[irow]
        fig.text(margin_left + 0.5*margin_x, 1. - margin_top - \
                0.5*margin_y - irow*(subplot_height + margin_y), label,\
                va='top', ha='left', fontsize=14,\
                bbox=dict(facecolor='white'))

    # Put some useful information on the title
    averaging_time = (times[-1] - times[0])/len(times)*navg/time_unit
    title = dirname_stripped + '     ' +\
            (r'$r/R_\odot\ =\ %0.3f$' %rval_to_plot) + '     ' +\
            (('t_avg = %.1f ' + time_label) %averaging_time)

    axs[0].set_title(title, **csfont)
    # Get ticks everywhere
    for ax in axs:
        plt.sca(ax)
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
