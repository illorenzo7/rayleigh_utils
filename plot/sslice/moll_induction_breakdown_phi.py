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
from plotcommon import axis_range
from sslice_util import plot_moll
from rayleigh_diagnostics import Shell_Slices
from get_sslice import smooth, prime
from varprops import texlabels

# Get command line arguments
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Data with Shell_Slices
radatadir = dirname + '/Shell_Slices/'

file_list, int_file_list, nfiles = get_file_lists(radatadir)

minmax = None
symlog = False
logscale = False
iiter = nfiles - 1 # by default plot the last iteration
ir = 0 # by default plot just below the surface
rval = None # can also find ir by finding the closest point
            # to a local radius divided by rsun
clon = 0
saveplot = False
ncol = 2
must_smooth = False

plotdir = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-ir':
        ir = int(args[i+1])
    elif arg == '-rval':
        rval = float(args[i+1])
    elif arg == '-smooth':
        dlon = int(args[i+1])
        print ("smoothing nonfield vars over %i degrees in lon." %dlon)
        must_smooth = True
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
    elif arg == '-clon':
        clon = float(args[i+1])
    elif arg == '-tag':
        rootname = args[i+1]
    elif arg == '-ncol':
        ncol = int(args[i+1])

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

if plotdir is None:
    plotdir = dirname + '/plots/'
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)

iter_val = int_file_list[iiter]
fname = file_list[iiter]

# Read in desired shell slice
a = Shell_Slices(radatadir + fname, '')

# compute induction breakdown

# fluid variables
vr = a.vals[:, :, :, a.lut[1], 0]
vt = a.vals[:, :, :, a.lut[2], 0]
vp = a.vals[:, :, :, a.lut[3], 0]

br = a.vals[:, :, :, a.lut[801], 0]
bt = a.vals[:, :, :, a.lut[802], 0]
bp = a.vals[:, :, :, a.lut[803], 0]

# fluid derivatives
# velocity field
dvrdr = a.vals[:, :, :, a.lut[10], 0]
dvtdr = a.vals[:, :, :, a.lut[11], 0]
dvpdr = a.vals[:, :, :, a.lut[12], 0]

dvrdt = a.vals[:, :, :, a.lut[37], 0]
dvtdt = a.vals[:, :, :, a.lut[38], 0]
dvpdt = a.vals[:, :, :, a.lut[39], 0]

dvrdp = a.vals[:, :, :, a.lut[46], 0]
dvtdp = a.vals[:, :, :, a.lut[47], 0]
dvpdp = a.vals[:, :, :, a.lut[48], 0]

# B field
dbrdr = a.vals[:, :, :, a.lut[810], 0]
dbtdr = a.vals[:, :, :, a.lut[811], 0]
dbpdr = a.vals[:, :, :, a.lut[812], 0]

dbrdt = a.vals[:, :, :, a.lut[837], 0]
dbtdt = a.vals[:, :, :, a.lut[838], 0]
dbpdt = a.vals[:, :, :, a.lut[839], 0]

dbrdp = a.vals[:, :, :, a.lut[846], 0]
dbtdp = a.vals[:, :, :, a.lut[847], 0]
dbpdp = a.vals[:, :, :, a.lut[848], 0]

# induction terms
sh = prime(a.vals[:, :, :, a.lut[1601], 0])
co = prime(a.vals[:, :, :, a.lut[1602], 0])
ad = prime(a.vals[:, :, :, a.lut[1603], 0])
ind_r = sh + co + ad

sh = prime(a.vals[:, :, :, a.lut[1606], 0])
co = prime(a.vals[:, :, :, a.lut[1607], 0])
ad = prime(a.vals[:, :, :, a.lut[1608], 0])
ind_t = sh + co + ad

sh = prime(a.vals[:, :, :, a.lut[1611], 0])
co = prime(a.vals[:, :, :, a.lut[1612], 0])
ad = prime(a.vals[:, :, :, a.lut[1613], 0])
ind_p = sh + co + ad

# will need some reference state stuff
rr = a.radius.reshape((1, 1, a.nr))
nt = a.ntheta
cost = a.costheta.reshape((1, nt, 1))
sint = a.sintheta.reshape((1, nt, 1))
cott = cost/sint

# now we need cost to be 1D
cost = a.costheta

# collect the terms
terms_r = [br, ind_r]
terms_r.append(prime(dvrdt*bt + vr*dbtdt))
terms_r.append(-prime(dvtdt*br + vt*dbrdt))
terms_r.append(-prime(dvpdp*br + vp*dbrdp))
terms_r.append(prime(dvrdp*bp + vr*dbpdp))
terms_r.append(prime(cott/rr*(vr*bt - vt*br)))

terms_t = [bt, ind_t]
terms_t.append(prime(dvtdp*bp + vt*dbpdp))
terms_t.append(-prime(dvpdp*bt + vp*dbtdp))
terms_t.append(-prime(dvrdr*bt + vr*dbtdr))
terms_t.append(prime(dvtdr*br + vt*dbrdr))
terms_t.append(prime(1./rr*(vt*br - vr*bt)))

terms_p = [bp, ind_p]
terms_p.append(prime(dvpdr*br + vp*dbrdr))
terms_p.append(prime(-vr*dbpdr - dvrdr*bp))
terms_p.append(prime(dvpdt*bt + vp*dbtdt))
terms_p.append(prime(-dvtdt*bp - vt*dbpdt))
terms_p.append(prime(1./rr*(vp*br - vr*bp)))

# term labels
labels_r = ['B_r', '[del X (v X B)]_r', '(d/dT)(vr*bt)',\
    '-(d/dT)(vt*br)', '-(d/dP)(vp*br)', '(d/dP)(vr*bp)',\
    '(cott/r)(vr*bt - vt*br)', 'sum last 5']

labels_t = ['B_theta', '[del X (v X B)]_theta', '(d/dP)(vt*bp)',\
    '-(d/dP)(vp*bt)', '-(d/dr)(vr*bt)', '(d/dr)(vt*br)',\
    '(1/r)(vt*br - vr*bt)', 'sum last 5']

labels_p = ['B_phi', '[del X (v X B)]_phi', '(d/dr)(vt*br)',\
    '-(d/dr)(vr*bt)', '(1/r)(d/dt)(vp*bt)', '-(1/r)(d/dt)(vt*bp)',\
    '(1/r)(vp*br - vr*bp)', 'sum last 5']

# Other info from the shell slice
t_loc = a.time[0]
# Find desired radius (by default ir=0--near outer surface)
if not rval is None:
    ir = np.argmin(np.abs(a.radius/rsun - rval))
rval = a.radius[ir]/rsun # in any case, this is the actual rvalue we get

# free some memory
del a

# figure dimensions
nplots = len(terms_r)
if ncol > nplots:
    ncol = nplots # (no point in having empty columns)
nrow = np.int(np.ceil(nplots/ncol))
fig_width_inches = 6.*ncol

# General parameters for main axis/color bar
margin_bottom_inches = 1./2.
margin_top_inches = 3./4. # margin for big title
margin_subplot_top_inches = 1/4 # margin to accommodate just subplot titles
margin_inches = 1./8.

subplot_width_inches = (fig_width_inches - (ncol + 1)*margin_inches)\
    /ncol
subplot_height_inches = 0.5*subplot_width_inches
fig_height_inches = margin_top_inches + nrow*(margin_bottom_inches +\
    subplot_height_inches + margin_subplot_top_inches)
fig_aspect = fig_height_inches/fig_width_inches

# "Non-dimensional" figure parameters
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_subplot_top = margin_subplot_top_inches/fig_height_inches

subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

# make (3) plots
for terms, labels, suf in [(terms_r, labels_r, '_r'), (terms_t, labels_t, '_theta'), (terms_p, labels_p, '_phi')]:
    rootname = 'induction_breakdown' + suf
    term_tot = np.zeros_like(terms[0])
    for i in range(2, nplots):
        term_tot += terms[i]
    terms.append(term_tot)

    print ("rms(exact - sum)/rms(exact) = ", rms(terms[1] - terms[-1])/rms(terms[1]))

    # Main figure
    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
        
    print('Plotting moll: ' + rootname + (', r/rsun = %0.3f (ir = %02i), '\
            %(rval, ir)) + 'iter ' + fname)

    # loop over desired vars and make plots
    for iplot in range(nplots):
        ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
        ax_bottom = 1. - margin_top - margin_subplot_top -\
                subplot_height - (iplot//ncol)*(subplot_height +\
                margin_subplot_top + margin_bottom)
        ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))

        if must_smooth and not labels[iplot] in ['B_r', 'B_theta', 'B_phi']:
            terms[iplot] = smooth(terms[iplot], dlon)
        field = terms[iplot][:, :, ir]
        # Display at terminal what we are plotting
        plot_moll(field, cost, fig=fig, ax=ax, minmax=minmax, clon=clon,\
                varname=rootname, symlog=symlog, logscale=logscale) 

        ax.set_title(labels[iplot], verticalalignment='bottom', **csfont)

        # Make title
        if iplot == 0:
            ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
            ax_delta_x = ax_xmax - ax_xmin
            ax_delta_y = ax_ymax - ax_ymin
            ax_center_x = ax_xmin + 0.5*ax_delta_x    

            if rotation:
                time_string = ('t = %.1f ' %(t_loc/time_unit)) + time_label +\
                        ' (1 ' + time_label + (' = %.2f days)'\
                        %(time_unit/86400.))
            else:
                time_string = ('t = %.3f ' %(t_loc/time_unit)) + time_label +\
                        ' (1 ' + time_label + (' = %.1f days)'\
                        %(time_unit/86400.))

            title = dirname_stripped +\
                '\n' + r'$\rm{Mollweide}$' + '     '  + time_string +\
                '\n' + (r'$r/R_\odot\ =\ %0.3f$' %rval)
            fig.text(ax_center_x, ax_ymax + 0.02*ax_delta_y +\
                    margin_subplot_top, title,\
                 verticalalignment='bottom', horizontalalignment='center',\
                 fontsize=10, **csfont)   

    plotdir = dirname + '/plots/moll/'
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)
    savefile = plotdir + 'moll_' + rootname + '_' + fname +\
            ('_rval%0.3f' %rval) + '.png'
    print ('Saving plot at ' + savefile)
    plt.savefig(savefile, dpi=300)