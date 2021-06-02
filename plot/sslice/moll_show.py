import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['rapl'])
from common import *
from plotcommon import *
from cla_util import *
from sslice_util import plot_moll
from rayleigh_diagnostics import Shell_Slices
#from get_sslice import get_sslice
from get_slice import get_slice, get_label
from varprops import texlabels

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Data with Shell_Slices
radatadir = dirname + '/Shell_Slices/'

file_list, int_file_list, nfiles = get_file_lists(radatadir, args)

# SPECIFIC ARGS for moll_show:
kwargs = dict({'iiter': nfiles - 1, 'irvals': np.array([0]), 'rvals': None, 'qvals': 'vr'})
# update these defaults from command-line
kwargs_moll = {**kwargs_contourf}
kwargs_moll['clon'] = 0.
del kwargs_moll['fig']
del kwargs_moll['ax']

kwargs = dotdict(update_kwargs(clas, kwargs))
kwargs_moll = dotdict(update_kwargs(clas, kwargs_moll))

iiter = kwargs.iiter
irval = kwargs.irvals[0] # these are going to be 1D arrays
rval = kwargs.rvals
qval = kwargs.qvals

# Get the baseline time unit
time_unit, time_label, rotation = get_time_unit(dirname)

# Read in desired shell slice
fname = file_list[iiter]
print ("reading sslice from " + radatadir + fname)
a = Shell_Slices(radatadir + fname, '')
if 'the_file' in clas:
    the_file = clas['the_file']
    print ("getting an averaged sslice from " + the_file)
    di = get_dict(the_file)
    a.vals = di['vals'][..., np.newaxis]
print ("done reading")

# figure dimensions
nplots = 1
sub_width_inches = 10.
sub_aspect = 1/2
margin_top_inches = 3/4 # larger top margin to make room for titles
margin_bottom_inches = 3/4
# larger bottom margin to make room for colorbar

# make plot
fig, axs, fpar = make_figure(nplots=nplots, sub_aspect=sub_aspect, margin_top_inches=margin_top_inches, margin_bottom_inches=margin_bottom_inches, sub_width_inches=sub_width_inches)
ax = axs[0, 0]

vals = get_slice(a, qval, dirname=dirname)

# Get local time (in seconds)
t_loc = a.time[0]

# Find desired radius (by default ir=0--near outer surface)
if not rval is None:
    rval = rval[0] # this was a 1D array
    irval = np.argmin(np.abs(a.radius/rsun - rval))
field = vals[:, :, irval]
rval = a.radius[irval]/rsun # in any case, this is the actual rvalue we get

# Display at terminal what we are plotting
print('Plotting moll: ' + qval + (', rval = %0.3f (irval = %02i), '\
        %(rval, irval)) + 'iter ' + fname)

plot_moll(field, a.costheta, fig=fig, ax=ax, **kwargs_moll)

# make title
if rotation:
    time_string = ('t = %.1f ' %(t_loc/time_unit)) + time_label
else:
    time_string = ('t = %.3f ' %(t_loc/time_unit)) + time_label
varlabel = texlabels[qval]

title = dirname_stripped +\
        '\n' + varlabel + '     '  + time_string +\
        '\n' + (r'$r/R_\odot\ =\ %0.3f$' %rval)
ax.set_title(title, va='bottom', **csfont)

# always show
plt.show()   
