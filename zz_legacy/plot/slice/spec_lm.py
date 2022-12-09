import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['rapl'])
from common import *
from plotcommon import *
from cla_util import *
from get_slice import get_label
from varprops import var_indices

# Get CLAs
args = sys.argv 
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# SPECIFIC ARGS for moll_show:
plot_specav_lm_kwargs_default = dict({'the_file': None, 'val_iter': 1e9, 'irvals': np.array([0]), 'rvals': None, 'varnames': np.array(['vr'])})
plot_specav_lm_kwargs_default.update(plot_spec_kwargs_default)
kw = update_dict(plot_specav_lm_kwargs_default, clas)
find_bad_keys(plot_specav_lm_kwargs_default, clas, 'plot_specav_lm', justwarn=True)
kw_plot_spec = update_dict(plot_spec_kwargs_default, clas)

# needs to be arrays
kw.irvals = make_array(kw.irvals)
kw.rvals = make_array(kw.rvals)
kw.varnames = make_array(kw.varnames)

# make plot directory if nonexistent
plotdir = my_mkdir(clas0['plotdir'] + 'specav_lm/')

# Read in spec data
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], 'Shell_Spectra')
print ('Reading Shell_Spectra data from '+ the_file)
di = get_dict(the_file)
vals = di['vals']
qv = di['qv']
lut = di['lut']

nell, nm, nr = np.shape(vals[:, :, :, 0])
lvals = np.arange(nell)
mvals = np.arange(nm)
rvals = di['rvals']

# Read in desired shell slice or average
# Data with Shell_Slices
if kw.the_file is None:
    kw.the_file = clas['the_file']
    print ("getting an averaged sslice from " + the_file)
    di = get_dict(the_file)
    a.vals = di['vals'][..., np.newaxis]
print ("done reading")

# get the rvals we want
if not kw.rvals is None: # irvals haven't been set directly
    if isall(kw.rvals):
        kw.irvals = np.arange(a.nr)
    else:
        kw.irvals = np.zeros_like(kw.rvals, dtype='int')
        for i in range(len(kw.rvals)):
            kw.irvals[i] = np.argmin(np.abs(a.radius/rsun - kw.rvals[i]))

# get the vars we want
if isall(kw.varnames):
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

    plt.sca(ax)
    # set bounds
#    plt.xlim(0.5, nell - 0.5)
#    plt.ylim(0.5, nm - 0.5)

    # label axes
    plt.xlabel(r'${\rm{spherical\ harmonic\ degree}}\ \ell$', fontsize=kw.fontsize)
    plt.ylabel(r'${\rm{azimuthal\ order}}\ m$', fontsize=kw.fontsize)

    # Get ticks everywhere
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    # Make title
    # Compute l_rms and m_rms
    l_rms = np.sum(field*lvals_2d)/np.sum(field)
    m_rms = np.sum(field*mvals_2d)/np.sum(field)

    # Make title
    title = dirname_stripped +\
        '\n' + r'$\rm{specav\_lm}$' + '     '  + time_string +\
        '\n' + varlabel + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %rval) +\
        '\n' + (r'$\ell_{\rm{rms}} = %.1f$' %l_rms) + '     ' +\
        (r'$m_{\rm{rms}} = %.1f$' %m_rms)
    fig.text(ax_center_x, ax_ymax + 0.02*ax_delta_y, title,\
         verticalalignment='bottom', horizontalalignment='center',\
         fontsize=fs, **csfont)   

    plt.savefig(plotdir + savename, dpi=300)
    if showplot:
        plt.show()
    plt.close()
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
