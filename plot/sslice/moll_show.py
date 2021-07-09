import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['rapl'])
from common import *
from plotcommon import *
from cla_util import *
from sslice_util import *
from rayleigh_diagnostics import Shell_Slices
from get_slice import get_slice, get_label

# Get CLAs
args = sys.argv 
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# SPECIFIC ARGS for moll_show:
moll_show_kwargs_default = dict({'val_iter': 1e9, 'irvals': np.array([0]), 'rvals': None, 'varnames': np.array(['vr'])})
moll_show_kwargs_default.update(plot_moll_kwargs_default)
kw = update_dict(moll_show_kwargs_default, clas)
find_bad_keys(moll_show_kwargs_default, clas, 'moll_show', justwarn=True)
kw_plot_moll = update_dict(plot_moll_kwargs_default, clas)

# needs to be arrays
kw.irvals = make_array(kw.irvals)
kw.rvals = make_array(kw.rvals)
kw.varnames = make_array(kw.varnames)

# make plot directory if nonexistent
plotdir = my_mkdir(clas0['plotdir'] + 'moll/')

# Read in desired shell slice or average
# Data with Shell_Slices
radatadir = dirname + '/Shell_Slices/'
fname = get_closest_file(radatadir, kw.val_iter)
print ("reading " + fname)
a = Shell_Slices(fname, '')
if 'the_file' in clas:
    the_file = clas['the_file']
    print ("getting an averaged sslice from " + the_file)
    di = get_dict(the_file)
    a.vals = di['vals'][..., np.newaxis]
print ("done reading")

# get the rvals we want
if not kw.rvals is None: # irvals haven't been set directly
    if kw.rvals == np.array(['all']):
        kw.irvals = np.arange(a.nr)
    else:
        kw.irvals = np.zeros_like(kw.rvals, dtype='int')
        for i in range(len(kw.rvals)):
            kw.irvals[i] = np.argmin(np.abs(a.radius/rsun - kw.rvals[i]))

# get the vars we want
if np.all(kw.varnames == 'all'): # remember varnames is an array now
    kw.varnames = get_default_varnames(dirname)

# plot dimensions
nplots = 1
sub_width_inches = 8
sub_aspect = 1/2
margin_top_inches = 3/8 # larger top margin to make room for titles
margin_bottom_inches = 1/2
# larger bottom margin to make room for colorbar

# loop over rvals/vars and make plots
print ("========================")
print ("about to plot:")
print ("irvals = ", kw.irvals)
print ("varnames = ", kw.varnames)
nfigures = len(kw.irvals)*len(kw.varnames)
print ("nfigures = ", nfigures)
print ("========================")

for varname in kw.varnames:
    # get the desired field variable
    vals = get_slice(a, varname, dirname=dirname)
    for irval in kw.irvals:
        field = vals[:, :, irval]
        rval = a.radius[irval]/rsun 

        # Display at terminal what we are plotting
        savename = 'moll_' + str(a.iters[0]).zfill(8) + '_' + varname + ('_rval%0.3f' %rval) + '.png' 

        # make plot
        fig, axs, fpar = make_figure(nplots=nplots, sub_width_inches=sub_width_inches, sub_aspect=sub_aspect, margin_top_inches=margin_top_inches, margin_bottom_inches=margin_bottom_inches)
        ax = axs[0, 0]
        plot_moll(field, a.costheta, fig, ax, **kw_plot_moll)

        # make title
        time_string = get_time_string(dirname, a.iters[0])
        varlabel = get_label(varname)

        title = dirname_stripped +\
                '\n' + varlabel + '     '  + time_string + '     ' +\
                (r'$r/R_\odot\ =\ %0.3f$' %rval) + '\n' +\
                ('clon = %4.0f' %kw.clon)
        ax.set_title(title, va='bottom', fontsize=default_titlesize)

        # save by default
        if clas0['saveplot']:
            print ("saving " + plotdir + savename)
            plt.savefig(plotdir + savename, dpi=300)
        # always show if nfigures is 1
        if nfigures == 1 and clas0['showplot']:
            print ("displaying " + savename[:-4])
            plt.show()   
        plt.close()
