# Routine to plot extended time-latitude cuts of B-field for the 
# "allthrees" sims
# Created: 11/10/2019
import matplotlib.pyplot as plt
from matplotlib import colors
SMALL_SIZE = 7
MEDIUM_SIZE = 9
BIGGER_SIZE = 11

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsiz e of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from common import get_file_lists, get_widest_range_file, strip_dirname,\
        rsun, get_dict, allthrees_start, get_exp, get_symlog_params,\
        sci_format
from plotcommon import axis_range
from rayleigh_diagnostics import GridInfo
from get_parameter import get_parameter

# Get the run directory on which to perform the analysis
dirname = sys.argv[1] 
dirname_stripped = strip_dirname(dirname)

# Data and plot directories
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# more defaults
minmax1 = None
minmax2 = None
dr = 0.03 # By default average energies over 10% the depth of the shell
        # (shell_depth ~ 0.27 Rsun)
gi = GridInfo('/altair/loma3853/dyn_nkeom3.0-1/grid_info', '')
rw = gi.rweights
saveplot = True
symlog = False
varname = 'bp' # by default, plot B_phi and KE_phi
nstd = None
linthresh = None

desired_rvals = [0.83] # by default, plot time-radius diagram for fields 
    # mid-CZ (units of solar radius)

# Get command-line arguments
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax1': # minmax1 is for time-latitude
        minmax1 = float(args[i+1]), float(args[i+2])
        print("Got CLA minmax1: ", minmax1)
    elif arg == '-minmax2':
        minmax2 = float(args[i+3]), float(args[i+4])
    elif arg == '-dr':
        dr = float(args[i+1])
    elif arg == '-rvals':
        string_desired_rvals = args[i+1].split()
        if string_desired_rvals == ['all']:
            desired_rvals = 'all'
        else:
            desired_rvals = []
            for j in range(len(string_desired_rvals)):
                desired_rvals.append(float(string_desired_rvals[j]))
    elif arg == '-var':
        varname = args[i+1]
    elif arg == '-symlog':
        symlog = True
    elif arg == '-nstd':
        nstd = float(args[i+1])
    elif arg == '-linthresh':
        linthresh = float(args[i+1])

# Get global rotation rate; this script fails for non-rotating models
angular_velocity = get_parameter(dirname, 'angular_velocity')
Prot = 2*np.pi/angular_velocity

the_file = get_widest_range_file(datadir, 'time-latitude')
print("Getting time-latitude data from ",  the_file)
di = get_dict(datadir + the_file)
vals = di['vals']

times = di['times']/Prot - allthrees_start
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

# Reshape the rweights to sum over time-radius Shell_Avgs slices
rw = rw.reshape((1, nr))

the_file = get_widest_range_file(datadir, 'trace_Shell_Avgs')
print("Getting shell avgs trace data from ",  the_file)
di_shtr = get_dict(datadir + the_file)
vals_shtr = di_shtr['vals']
lut_shtr = di_shtr['lut']
times_shtr = di_shtr['times']/Prot - allthrees_start

# Parameters to divide the extended trace into rows
t_begin = times[0]
t_end = times[-1]
Delta_t = t_end - t_begin
interval = 2220.26
nrows = int(np.ceil(Delta_t/interval))

# Make figure, stretching time-latitude plots over multiple rows, with enery trace beneath
fig_width_inches = 7+1/4
margin_inches = 1/8
margin_vert_inches = 1/4 # space for time label
margin_bottom_inches = 3/8 # space for x-axis label
if symlog:
    margin_top_inches = 5/8 # space for colorbar
else:
    margin_top_inches = 1/8
margin_left_inches = 5/8 # space for latitude label
subplot_width_inches = fig_width_inches - margin_inches - margin_left_inches
subplot_pad_inches = 0. # space between time-latude and energy-trace pairs
two_subplots_height_inches = 1.3
tl_height_inches = 0.9
etr_height_inches = two_subplots_height_inches - tl_height_inches - subplot_pad_inches

fig_height_inches = nrows*two_subplots_height_inches +\
        (nrows - 1)*margin_vert_inches + margin_bottom_inches +\
        margin_top_inches
fig_aspect = fig_height_inches/fig_width_inches

margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_vert = margin_vert_inches/fig_height_inches
tl_height = tl_height_inches/fig_height_inches
etr_height = etr_height_inches/fig_height_inches
subplot_width = subplot_width_inches/fig_width_inches
margin_left = margin_left_inches/fig_width_inches
margin_bottom = margin_bottom_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches

two_subplots_height = two_subplots_height_inches/fig_height_inches

# Find the r-values we want
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

# Get raw traces of variable
q_dict_tl = {'br': 801, 'bt': 802, 'bp': 803}
q_dict_shtr = {'br': 1102, 'bt': 1103, 'bp': 1104}
alabels = {'br': r'$\langle B_r\rangle$', 'bt': r'$\langle B_\theta\rangle$',\
        'bp': r'$\langle B_\phi\rangle$'}
blabels = {'br': r'$\rm{ME_r}$', 'bt': r'$\rm{ME_\theta}$', \
        'bp': r'$\rm{ME_\phi}$'}

iq = np.argmin(np.abs(qvals - q_dict_tl[varname]))
quant = vals[:, :, :, iq]
quant_energy = vals_shtr[:, :, lut_shtr[q_dict_shtr[varname]]]

# Loop over the desired radii and save plots
for i in range(len(i_desiredrvals)):
    i_desiredrval = i_desiredrvals[i]
    print("=====================================")
    print ("Plotting ", varname, ", ir = ", i_desiredrval)
    rval_to_plot = rvals_to_plot[i]
    # Time-latitude of quantity
    quant_loc = quant[:, :, i_desiredrval]
    print("minmax1 is ", minmax1)
    print("len(i_desired_rvals) is ", len(i_desiredrvals))
    if minmax1 is None or len(i_desiredrvals) > 1:
        if nstd is None:
            if symlog:
                nstd = 9.
            else:
                nstd= 3.
        minmax1 = -nstd*np.std(quant_loc), nstd*np.std(quant_loc)
        print ("Setting time-latitude minmax to (%1.2e, %1.2e), (nstd = %.1f)"\
                %(minmax1[0], minmax1[1], nstd))
    # Associated energy, in range r_loc \pm dr/2
    ir1 = np.argmin(np.abs(rr/rsun - (rval_to_plot + dr/2.)))
    ir2 = np.argmin(np.abs(rr/rsun - (rval_to_plot - dr/2.)))
    print("Averaging energy over (ir1, ir2) = ", ir1, ir2)
    print ("or (rmin/rsun, rmax/rsun) = (%1.3f, %1.3f)" %(rr[ir2]/rsun,\
            rr[ir1]/rsun))
    quant_energy_loc = np.sum((quant_energy*rw)[:, ir1:ir2+1], axis=1)/\
            np.sum(rw[0, ir1:ir2+1])
    quant_energy_exp = get_exp(np.max(quant_energy_loc))
    quant_energy_loc /= 10.0**quant_energy_exp
    if minmax2 is None or len(i_desiredrvals) > 1:
        minmax2 = 0., 1.05*np.max(quant_energy_loc)
        print ("Setting energy minmax to (%1.2e, %1.2e)"\
                %(minmax2[0], minmax2[1]*quant_energy_exp))

    if symlog:
        # get symlog parameters
        linthresh_tmp, linscale = get_symlog_params(quant_loc, field_max=minmax1[1])
        if linthresh is None:
            linthresh = linthresh_tmp

    # Make the figure, spanning multiple rows
    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    axs = []
    for irow in range(nrows): 
        axs.append(fig.add_axes((margin_left, 1 - margin_top - tl_height -\
                                 irow*(margin_vert + two_subplots_height),\
                    subplot_width, tl_height)))                 

        axs.append(fig.add_axes((margin_left,\
                1 - margin_top - two_subplots_height -\
                irow*(margin_vert + two_subplots_height),\
                subplot_width, etr_height)))

        desired_t1 = t_begin + interval*irow
        desired_t2 = t_begin + interval*(irow + 1)
        it1 = np.argmin(np.abs(times - desired_t1))
        it2 = np.argmin(np.abs(times - desired_t2))
        
        it1_shtr = np.argmin(np.abs(times_shtr - desired_t1))
        it2_shtr = np.argmin(np.abs(times_shtr - desired_t2))
        
        times_interval = times[it1:it2+1]
        times_interval_shtr = times_shtr[it1_shtr:it2_shtr+1]
        
        times_2d, tt_lat_2d = np.meshgrid(times_interval, tt_lat, indexing='ij')

        quant_interval = quant_loc[it1:it2+1, :]
        quant_energy_interval = quant_energy_loc[it1_shtr:it2_shtr + 1]

        if symlog:
            norm = colors.SymLogNorm(linthresh=linthresh,\
                linscale=linscale, vmin=minmax1[0], vmax=minmax1[1])
            im = axs[2*irow].pcolormesh(times_2d, tt_lat_2d, quant_interval,\
                    cmap='RdYlBu_r', norm=norm)
        else:
            im = axs[2*irow].pcolormesh(times_2d, tt_lat_2d, quant_interval,\
                    cmap='RdYlBu_r', vmin=minmax1[0], vmax=minmax1[1])

        axs[2*irow + 1].plot(times_interval_shtr, quant_energy_interval, 'k',\
                linewidth=0.5)

        # Set reasonable latitude labels
        axs[2*irow].set_yticks([-60, -30, 0, 30, 60])
#        axs[2*irow].set_yticklabels([])
#        axs[2*irow + 1].set_yticklabels([])
    #    axs[2*irow + 1].set_yticks([-60, -30, 0, 30, 60])   

        # Get ticks everywhere
        plt.sca(axs[2*irow])
        plt.minorticks_on()
        plt.tick_params(top=True, right=True, direction='in', which='both')

        plt.sca(axs[2*irow + 1])
        plt.minorticks_on()
        plt.tick_params(top=True, right=True, direction='in', which='both')

        # Get rid of tick labels on "top subplot" panels
        axs[2*irow].set_xticklabels([])

        # Label y-axis
        label_fs = 7
        axs[2*irow].set_ylabel(r'$\rm{latitude\ (^\circ)}$',\
                fontsize=label_fs)
#        axs[2*irow + 1].set_ylabel(r'$\rm{[10^{%i}\ erg\ cm^{-3}]}$'\
        axs[2*irow + 1].set_ylabel(r'$\rm{[\times10^{%i}\ cgs]}$'\
                %quant_energy_exp, fontsize=label_fs, labelpad=15)

        axs[2*irow].set_xlim(desired_t1, desired_t2)
        axs[2*irow].set_ylim(-90., 90.)
       
        axs[2*irow + 1].set_xlim(desired_t1, desired_t2)
        axs[2*irow + 1].set_ylim(minmax2[0], minmax2[1])
        
    #    axs[2*irow + 1].ticklabel_format(scilimits = (-3,4), useMathText=True)

    # Label with saturation values
    fs = 9
    axs[0].text(100., 75., '(a) ' + alabels[varname],\
            va='top', ha='left', fontsize=fs)
    axs[0].text(2135., 75., r'$\pm %1.2f\ kG$' %(minmax1[1]/1000),\
            va='top', ha='right', fontsize=fs)
    axs[0].text(100., -75., r'$r/R_\odot = %0.3f$'\
            %rvals_sampled[i_desiredrval], va='bottom', ha='left',\
            fontsize=fs)

    axs[1].text(100, 0.8*minmax2[1], '(b) ' + blabels[varname],\
            va='top', ha='left', fontsize=fs)

    # Label the axes
    axs[2*nrows - 1].set_xlabel(r'$t\ (P_{\rm{rot}})$')

    # Label with the case name in upper left
    if dirname_stripped == 'dyn_nkeom3.0-alldata':
        label = 'Case D3-1'
    elif dirname_stripped == 'dyn_nk3.0_e1.5_om3.0-alldata':
        label = 'Case D3-2'
    elif dirname_stripped == 'dyn_nk3.0_e0.8_om3.0-alldata':
        label = 'Case D3-4'

    if i_desiredrval == 4 and varname == 'bp' and \
            dirname_stripped == 'dyn_nkeom3.0-alldata':
        # Draw some arrows where asymmetric cycle is "born" from symmetric cycle
        axs[0].arrow(2125, -85, 0, 40, head_width=10, head_length=10,\
                fc='k', ec='k')
        axs[4].arrow(4730, -85, 0, 40, head_width=10, head_length=10,\
                fc='k', ec='k')
        axs[6].arrow(8400, -85, 0, 40, head_width=10, head_length=10,\
                fc='k', ec='k')
#        axs[8].arrow(9265, 85, 0, -40, head_width=10, head_length=10,\
#                fc='k', ec='k')

        # Mark region in blow-up with vertical dashed lines
        lw = 0.7
        axs[0].plot(1100 + np.zeros(100), np.linspace(-90, 90, 100),\
                'k--', linewidth=lw)
        axs[0].plot(1300 + np.zeros(100), np.linspace(-90, 90, 100),\
                'k--', linewidth=lw)

    if i_desiredrval == 4 and varname == 'bp' and \
            dirname_stripped == 'dyn_nkeom3.0-alldata':
         # Mark region in blow-up with vertical dashed lines
        lw = 0.7
        axs[2].plot(3150 + np.zeros(100), np.linspace(-90, 90, 100),\
                'k--', linewidth=lw)
        axs[2].plot(3450 + np.zeros(100), np.linspace(-90, 90, 100),\
                'k--', linewidth=lw)
      
    # Set up the colorbar
    if symlog:
        cbar_width = subplot_width/3
        cbar_aspect = 1/20
        cbar_height = cbar_width*cbar_aspect/fig_aspect
        cbar_left = margin_left + 0.5*(subplot_width - cbar_width)
        cbar_bottom = 1 - margin_y - cbar_height
        cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))

        cbar = plt.colorbar(im, cax=cax, orientation='horizontal')

        nlin = 5
        nlog = 6
        lin_ticks = np.linspace(-linthresh, linthresh, nlin)
        log_ticks1 = np.linspace(minmax1[0], -linthresh, nlog, endpoint=False)
        log_ticks2 = -log_ticks1[::-1]
        ticks = np.hstack((log_ticks1, lin_ticks, log_ticks2))
        nticks = nlin + 2*nlog
        cbar.set_ticks(ticks)
        #ticklabels = np.zeros(nticks, dtype='str')
        ticklabels = []
        for i in range(nticks):
            ticklabels.append(r'')
        #ticklabels = np.array(ticklabels)    
        ticklabels[0] = sci_format(minmax1[0])
        ticklabels[nlog] = sci_format(-linthresh)
        ticklabels[nticks//2] = r'$0$'
        ticklabels[nlog + nlin - 1] = sci_format(linthresh)
        ticklabels[nticks - 1] = sci_format(minmax1[1])
        cbar.set_ticklabels(ticklabels)
#            cbticks = [minmax1[0], -linthresh, 0, linthresh, minmax1[1]]
    #        ticks = minmax1[1]*(2.0*norm(cbticks) - 1.0)

    #        cbar.set_ticks(ticks)

#            cbar.set_ticks(cbticks)
#            cbar.set_ticklabels([sci_format(minmax1[0]), sci_format(-linthresh),\
#                    r'$0$', sci_format(linthresh), sci_format(minmax1[1])])

        fig.text(cbar_left + cbar_width + 1/8/fig_width_inches,\
            cbar_bottom + 0.5*cbar_height, r'$\rm{G}$', va='center',\
            **csfont, fontsize=7)

        # Put case label in upper left
        fig.text(margin_left, 1 - margin_top + 1/16/fig_height_inches,\
                label, va='bottom', ha='left', **csfont, fontsize=fs)

    # Save it in the "tl[/symlog]" subdirectory of the "figure_set" directory

    savedir = '/home5/loma3853/Desktop/Publications/allthrees/figure_set/tl/'
    if symlog:
        savedir += 'symlog/'
    savename = dirname_stripped + '_whole_' + varname + ('_trace_rval%0.3f.png'

            %rvals_sampled[i_desiredrval])

    print ("Saving figure at ", savedir + savename)
    print("=====================================")
    plt.savefig(savedir + savename, dpi=300)
    plt.close()
