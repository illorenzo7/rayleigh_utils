import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
from common import *
from plotcommon import *
from cla_util import *
from slice_util import spec_2D_fig_dimensions, plot_spec_2D

# Get CLAs
args = sys.argv 
clas0, clas = read_clas(args)
dirname = clas0.dirname
dirname_stripped = strip_dirname(dirname)

# SPECIFIC ARGS
kwargs_default = dotdict(dict({'the_file': None, 'irvals': np.array([0]), 'rvals': None, 'varnames': np.array(['vr'])}, 'mode': 'lpower', 'lval': None, 'mval': None))
# "mode" can be: lpower, mpower, l, m

# many kwargs
kwargs_default.update(plot_spec_2D_kwargs_default)
make_figure_kwargs_default.update(spec_2D_fig_dimensions)
kwargs_default.update(make_figure_kwargs_default)
find_bad_keys(kwargs_default, clas, 'plot/timetrace/tspec', justwarn=True)

kw = update_dict(kwargs_default, clas)
kw_plot_spec_2D = update_dict(plotting_func_kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)

# get the rvals we want
irvals = kw.irvals
if not kw.rvals is None: # irvals haven't been set directly
    if np.all(kw.rvals == 'all'):
        irvals = np.arange(a0.nr)
    else:
        irvals = np.zeros_like(kw.rvals, dtype='int')
        for i in range(len(kw.rvals)):
            irvals[i] = np.argmin(np.abs(a0.radius/rsun - kw.rvals[i]))

# and the qvals
qvals = make_array(kw.qvals)

# everything must be array
irvals = make_array(irvals)

# needs to be arrays
kw.irvals = make_array(kw.irvals)
kw.rvals = make_array(kw.rvals)
kw.varnames = make_array(kw.varnames)

# get data
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], dataname)

datatype = 'timelon_clat' + lat_format(clat) + '_dlat%03.0f' %dlat

for qval in 
    datadir = dirname + '/data/'
    if kw.the_file is None:
        kw.the_file = get_widest_range_file(datadir, dataname)
        iter1, iter2 = get_iters_from_file(kw.the_file)
        print ("plotting average file: " + kw.the_file)
        di = get_dict(kw.the_file)
        a0.vals = di['vals'][..., np.newaxis]

# get the rvals we want
if not kw.rvals is None: # irvals haven't been set directly
    if np.all(kw.rvals == 'all'):
        kw.irvals = np.arange(a0.nr)
    else:
        kw.irvals = np.zeros_like(kw.rvals, dtype='int')
        for i in range(len(kw.rvals)):
            kw.irvals[i] = np.argmin(np.abs(a0.radius/rsun - kw.rvals[i]))

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
        print ("got here")
        print ("maxabs vr = ", np.max(np.abs(a.vals[:,:,0,a.lut[1],0])))
    else:
        a = reading_func(radatadir + fname, '')
    for varname in kw.varnames:
        # get the desired field variable
        vals = get_slice(a, varname, dirname=dirname)
        if plottype == 'speclm' and not kw.av:
            vals = np.abs(vals)**2

        for irval in kw.irvals:
            field = vals[:, :, irval]
            rval = a.radius[irval]/rsun 

            # Display at terminal what we are plotting
            if kw.av:
                savename = basename + '_' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8)
            else:
                savename = basename + '_' + str(a.iters[0]).zfill(8)
            savename += ('_' + varname + ('_rval%0.3f' %rval) + '.png')

            # make plot
            fig, axs, fpar = make_figure(**kw_make_figure)
            ax = axs[0, 0]
            if plottype == 'moll':
                plotting_args = field, a.costheta, fig, ax
            if plottype == 'speclm':
                plotting_args = field, fig, ax
                kw_plotting_func.cbar_pos = 'right'

            plotting_func(*plotting_args, **kw_plotting_func)

            # make lebels
            ax.set_xlabel('sph. harm. deg. (l)', fontsize=default_labelsize)
            ax.set_ylabel('azimuthal wavenumber (m)', fontsize=default_labelsize)

            # make title
            if kw.av:
                time_string = get_time_string(dirname, iter1, iter2, oneline=True)
            else:
                time_string = get_time_string(dirname, a.iters[0])
            varlabel = get_label(varname)

            if plottype == 'moll':
                slice_info = varlabel + 5*' ' + (r'$r/R_\odot\ =\ %0.3f$' %rval) + 5*' ' + ('clon = %4.0f' %kw.clon)
            if plottype == 'speclm':
                slice_info = varlabel + 5*' ' + (r'$r/R_\odot\ =\ %0.3f$' %rval)

            title = dirname_stripped + '\n' + slice_info + '\n' + time_string
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
print (buff_line)
