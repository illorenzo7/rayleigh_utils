import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
sys.path.append(os.environ['rapp'])
from common import *
from plotcommon import *
from sslice_util import plot_moll
from rayleigh_diagnostics import Shell_Slices
from get_slice import get_slice, get_label
from cla_util import *

# Get CLAs
args = sys.argv 
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# SPECIFIC ARGS for moll_show:
kwargs = dict({'val_iter': 1e9, 'irvals': np.array([0]), 'rvals': None, 'varname': 'vr'})
kwargs_moll = {**kwargs_contourf}
kwargs_moll['clon'] = 0.
kwargs_moll = dotdict(update_kwargs(clas, kwargs_moll))

# update these defaults from command-line
if 'di_trans' in clas:
    clas['val_iter'] = clas['di_trans']['val_iter']
kwargs = dotdict(update_kwargs(clas, kwargs))
irval = kwargs.irvals[0] # these are going to be 1D arrays
rvals = kwargs.rvals
varname = kwargs.varname

# make plot directory if nonexistent
plotdir = my_mkdir(clas0['plotdir'] + 'moll/rvals/')

# Read in desired shell slice or average
# Data with Shell_Slices
radatadir = dirname + '/Shell_Slices/'
fname = get_closest_file(radatadir, kwargs.val_iter)
print ("reading " + fname)
a = Shell_Slices(fname, '')
if 'the_file' in clas:
    the_file = clas['the_file']
    print ("getting an averaged sslice from " + the_file)
    di = get_dict(the_file)
    a.vals = di['vals'][..., np.newaxis]
print ("done reading")

# get the desired field variable
vals = get_slice(a, varname, dirname=dirname)

# plot dimensions
nplots = 1
sub_width_inches = 6
sub_aspect = 1/2
margin_top_inches = 3/8 # larger top margin to make room for titles
margin_bottom_inches = 1/2
# larger bottom margin to make room for colorbar

# Loop over rvals and make plots
if rvals is None:
    irvals = np.arange(a.nr) # loop over everything
else:
    irvals = np.zeros_like(rvals, dtype='int')
    for i in range(len(rvals)):
        irvals[i] = np.argmin(np.abs(a.radius/rsun - rvals[i]))

for ir in irvals:
    rval = a.radius[ir]/rsun
    field = vals[:, :, ir]

    savename = 'moll_' + varname + '_iter' + str(a.iters[0]).zfill(8) +\
            ('_rval%0.3f' %rval) + '.png'
    # Display at terminal what we are plotting
    # Display at terminal what we are plotting
    print('Plotting moll: ' + varname + (', rval = %0.3f (irval = %02i), ' %(rval, irval)) + 'iter %08i' %a.iters[0])

    # make plot
    fig, axs, fpar = make_figure(nplots=nplots, sub_width_inches=sub_width_inches, sub_aspect=sub_aspect, margin_top_inches=margin_top_inches, margin_bottom_inches=margin_bottom_inches)
    ax = axs[0, 0]
    plot_moll(field, a.costheta, fig, ax, **kwargs_moll)

    # make title
    time_string = get_time_string(dirname, a.iters[0])
    varlabel = get_label(varname)

    title = dirname_stripped +\
            '\n' + varlabel + '     '  + time_string + '     ' +\
            (r'$r/R_\odot\ =\ %0.3f$' %rval) + '\n' +\
            ('clon = %4.0f' %kwargs_moll['clon'])
    ax.set_title(title, va='bottom', fontsize=default_titlesize)

    plt.savefig(plotdir + savename, dpi=300)
    plt.close()
