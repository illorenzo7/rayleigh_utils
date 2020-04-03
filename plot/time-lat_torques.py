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
from time_scales import compute_Prot, compute_tdt

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]

# Data and plot directories
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
nosave = False
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)
dirname_stripped = strip_dirname(dirname)

# Find the time/latitude file(s) the data directory. If there are 
# multiple, by default choose the one with widest range in the trace.
the_file = get_widest_range_file(datadir, 'time-latitude_torques')

# more defaults
minmax = None
xminmax = None
tag = ''

desired_rvals = [0.692] # by default, plot time-radius diagram for fields 
    # mid-RZ (units of solar radius)
navg = 11 # by default average over 11 AZ_Avgs files for each time

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
    elif arg == '-tag':
        tag = '_' + args[i+1]

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

times_trace = times[over2:niter - over2]/time_unit

# Make meshgrid of time/radius
# Take into account if user specified xmin, xmax
if xminmax is None: # By default use all times available
    it1 = 0
    it2 = niter - navg
else:
    it1 = np.argmin(np.abs(times_trace - xmin))
    it2 = np.argmin(np.abs(times_trace - xmax))
times_trace = times_trace[it1:it2+1]
bp_all = bp_all[it1:it2+1]
trans_tot_all = trans_tot_all[it1:it2+1]
trans_mm_all = trans_mm_all[it1:it2+1]
diff_all = diff_all[it1:it2+1]
times2, tt_lat2 = np.meshgrid(times_trace, tt_lat, indexing='ij')

# Loop over the desired radii and save plots
for i in range(len(i_desiredrvals)):
    i_desiredrval = i_desiredrvals[i]
    rval_to_plot = rvals_to_plot[i]
    bp_trace = bp_all[:, :, i_desiredrval]
    trans_tot_trace = trans_tot_all[:, :, i_desiredrval]
    trans_mm_trace = trans_mm_all[:, :, i_desiredrval]
    trans_pp_trace = trans_tot_trace - trans_mm_trace
    diff_trace = diff_all[:, :, i_desiredrval]
    tot_trace = trans_tot_trace + diff_trace
  
    # Make appropriate file name to save
    savename = dirname_stripped + '_time-lat_induction_p_' +\
        ('rval%0.3f_' %rval_to_plot) + str(iter1).zfill(8) + '_' +\
        str(iter2).zfill(8) + tag + '.png'

    if not user_specified_minmax:
        std_bp = np.std(bp_trace)
        my_min_bp, my_max_bp = -3.*std_bp, 3.*std_bp

        std_diff = np.std(diff_trace)
        std_trans_tot = np.std(trans_tot_trace)
        std_trans_mm = np.std(trans_mm_trace)
        std_trans_pp = np.std(trans_pp_trace)
        std_tot = np.std(tot_trace)
        std_max = max(std_diff, std_trans_tot, std_trans_mm, std_trans_pp,\
                std_tot)
        my_min_induction, my_max_induction = -3.*std_max, 3.*std_max
     
    # Create figure with  5 panels in a row (phi induction terms:
    # mean-mean transport, prime-prime transport, total transport, diffusion,
    # actual bphi )
    fig, axs = plt.subplots(5, 1, figsize=(12, 8), sharex=True, sharey=True)
    ax1 = axs[0]; ax2 = axs[1]; ax3 = axs[2]; ax4 = axs[3]; ax5 = axs[4]

    im1 = ax1.pcolormesh(times2, tt_lat2, trans_mm_trace,\
            vmin=my_min_induction, vmax=my_max_induction, cmap='RdYlBu_r')
    im2 = ax2.pcolormesh(times2, tt_lat2, trans_pp_trace,\
            vmin=my_min_induction, vmax=my_max_induction, cmap='RdYlBu_r')
    im3 = ax3.pcolormesh(times2, tt_lat2, diff_trace,\
            vmin=my_min_induction, vmax=my_max_induction, cmap='RdYlBu_r')
    im4 = ax4.pcolormesh(times2, tt_lat2, tot_trace,\
            vmin=my_min_induction, vmax=my_max_induction, cmap='RdYlBu_r')
    im5 = ax5.pcolormesh(times2, tt_lat2, bp_trace,\
            vmin=my_min_bp, vmax=my_max_bp, cmap='RdYlBu_r')

    # Put colorbar next to all plots (possibly normalized separately)
    # First make room and then find location of subplots
    plt.subplots_adjust(left=0.1, right=0.85, wspace=0.03, top=0.9)

    # First, trans_mm:
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

    # Next, db dt:
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

    # Next, B_phi:
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
    label1 = r'$[\nabla\times(\overline{\mathbf{v}}\times\overline{\mathbf{B}})]_\phi$'
    ax1.text(ax_xmin + 1.01*ax_Dx, 0., label1, va='center',\
            rotation=270)
    label1 = r'$[\nabla\times(\overline{\mathbf{v}^\prime\times\mathbf{B}^\prime})]_\phi$'
    ax2.text(ax_xmin + 1.01*ax_Dx, 0.,  label1, va='center',\
            rotation=270)
    label1 = r'$[-\nabla\times(\eta(r)\nabla\times\overline{\mathbf{B}})]_\phi$'
    ax3.text(ax_xmin + 1.01*ax_Dx, 0., label1, va='center',\
            rotation=270)
    label1 = r'$\partial B_\phi/\partial t$'
    ax4.text(ax_xmin + 1.01*ax_Dx, 0.,  label1, va='center',\
            rotation=270)
    ax5.text(ax_xmin + 1.01*ax_Dx, 0.,  r'$\overline{B_\phi}$', va='center',\
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
