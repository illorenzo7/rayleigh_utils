import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['rapl'])
from common import *
from plotcommon import *
from cla_util import *
from slice_util import *
from rayleigh_diagnostics import Shell_Slices, Shell_Spectra
#from get_slice import get_slice, get_label

# Get CLAs
args = sys.argv 
clas0, clas = read_clas(args)
dirname = clas0.dirname
dirname_stripped = strip_dirname(dirname)

# default fig dimensions
moll_fig_dimensions = dict({'sub_width_inches': 6, 'sub_aspect': 1/2, 'sub_margin_left_inches': default_margin, 'sub_margin_top_inches': 1/2, 'sub_margin_bottom_inches': 1/2, 'margin_top_inches': 1/4})

spec_lm_fig_dimensions = dict({'sub_width_inches': 6, 'sub_aspect': 1, 'sub_margin_left_inches': default_margin_ylabel, 'sub_margin_top_inches': 1/2, 'sub_margin_bottom_inches': 1/2, 'sub_margin_right_inches': 7/8, 'margin_top_inches': 1/4})

# SPECIFIC ARGS
kwargs_default = dotdict(dict({'the_file': None, 'plottype': 'moll', 'av': False, 'val_iter': int(1e9), 'irvals': np.array([0]), 'rvals': None, 'varnames': np.array(['vr'])}))
# this guy need to update right away to choose fig dimensions
if 'plottype' in clas:
    kwargs_default.plottype = clas.plottype
   
if kwargs_default.plottype == 'moll':
    fig_dimensions = moll_fig_dimensions
    plotting_func = plot_moll 
    plotting_func_kwargs_default = plot_moll_kwargs_default
    dataname = 'Shell_Slices'
    reading_func = Shell_Slices
if kwargs_default.plottype == 'speclm':
    fig_dimensions = spec_fig_dimensions
    plotting_func = plot_spec_lm
    plotting_func_kwargs_default = plot_spec_lm_kwargs_default
    dataname = 'Shell_Spectra'
    reading_func = Shell_Spectra

kwargs_default.update(plotting_func_kwargs_default)
make_figure_kwargs_default.update(fig_dimensions)
kwargs_default.update(make_figure_kwargs_default)

find_bad_keys(kwargs_default, clas, 'plot/slice/generic', justwarn=True)

kw = update_dict(kwargs_default, clas)
kw_plotting_func = update_dict(plotting_func_kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)

# needs to be arrays
kw.irvals = make_array(kw.irvals)
kw.rvals = make_array(kw.rvals)
kw.varnames = make_array(kw.varnames)

# make plot directory if nonexistent
basename = kw.plottype
if kw.av:
    basename += 'av'
plotdir = my_mkdir(clas0['plotdir'] + basename + '/')

# get shell slice, range (for movie), or average
onefile = True
for arg in args:
    if arg in range_options:
        onefile = False
if onefile:
    args += ['--iter', str(kw.val_iter)]

# Get desired file names in datadir and their integer counterparts
radatadir = dirname + '/' + dataname + '/'
file_list, int_file_list, nfiles = get_file_lists(radatadir, args)

# need one of these no matter what
print ("reading " + dataname + '/' + file_list[0])
a0 = reading_func(radatadir + file_list[0], '')
print ("done reading")

if kw.av: # use time-averaged file
    datadir = dirname + '/data/'
    if kw.the_file is None:
        kw.the_file = get_widest_range_file(datadir, dataname)
        iter1, iter2 = get_iters_from_file(kw.the_file)
        print ("plotting average file: " + kw.the_file)
        di = get_dict(kw.the_file)
        a0.vals = di['vals'][..., np.newaxis]
else:
    print ('plotting %i %s files: %s through %s'\
        %(nfiles, dataname, file_list[0], file_list[-1]))

# get the rvals we want
if not kw.rvals is None: # irvals haven't been set directly
    if np.all(kw.rvals == 'all'):
        kw.irvals = np.arange(a.nr)
    else:
        kw.irvals = np.zeros_like(kw.rvals, dtype='int')
        for i in range(len(kw.rvals)):
            kw.irvals[i] = np.argmin(np.abs(a.radius/rsun - kw.rvals[i]))

# get the vars we want
if np.all(kw.varnames == 'all'): # remember varnames is an array now
    kw.varnames = get_default_varnames(dirname)

# loop over rvals/vars and make plots
print (buff_line)
print ("about to plot:")
print ("irvals = ", kw.irvals)
print ("varnames = ", kw.varnames)
print ("nfiles = ", nfiles)
nfigures = len(kw.irvals)*len(kw.varnames)*nfiles
print ("nfigures = ", nfigures)
print (buff_line)

for fname in file_list:
    if fname == file_list[0]:
        a = a0
    else:
        a = reading_func(radatadir + fname, '')
    for varname in kw.varnames:
        # get the desired field variable
        vals = get_slice(a, varname, dirname=dirname)
        for irval in kw.irvals:
            field = vals[:, :, irval]
            rval = a.radius[irval]/rsun 

            # Display at terminal what we are plotting
            if kw.av:
                savename = basename + '_' + str(iter1).zfill(8) + '_' + str(iter1).zfill(8)
            else:
                savename = basename + str(a.iters[0]).zfill(8)
            savename += ('_' + varname + ('_rval%0.3f' %rval) + '.png')

            # make plot
            fig, axs, fpar = make_figure(**kw_make_figure)
            ax = axs[0, 0]
            if kw.plottype == 'moll':
                plotting_args = field, a.costheta, fig, ax
            if kw.plottype == 'speclm':
                plotting_args = field, fig, ax
            plotting_func(*plotting_args, **kw_plotting_func)

            # make title
            if kw.av:
                time_string = get_time_string(dirname, iter1, iter2, oneline=True)
            else:
                time_string = get_time_string(dirname, a.iters[0])
            varlabel = get_label(varname)

            title = dirname_stripped + '\n' +\
                varlabel + 5*' ' + (r'$r/R_\odot\ =\ %0.3f$' %rval) + 5*' ' + ('clon = %4.0f' %kw.clon) + '\n' +\
                time_string
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
