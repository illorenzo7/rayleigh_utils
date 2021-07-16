import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
sys.path.append(os.environ['rapl'] + '/slice')
from common import *
from plotcommon import *
from cla_util import *
from slice_util import spec_2D_fig_dimensions, plot_spec_2D, plot_spec_2D_kwargs_default

# Get CLAs
args = sys.argv 
clas0, clas = read_clas(args)
dirname = clas0.dirname
dirname_stripped = strip_dirname(dirname)

# SPECIFIC ARGS
kwargs_default = dotdict(dict({'the_file': None, 'irvals': np.array([0]), 'rvals': None, 'qvals': np.array([1]), 'modes': np.array(['lpower']), 'lval': None, 'mval': None}))
# "modes" can be: lpower, mpower, l, m, sect, or combinations thereof

# many kwargs
kwargs_default.update(plot_spec_2D_kwargs_default)
make_figure_kwargs_default.update(spec_2D_fig_dimensions)
kwargs_default.update(make_figure_kwargs_default)
find_bad_keys(kwargs_default, clas, 'plot/timetrace/tspec', justwarn=True)

kw = update_dict(kwargs_default, clas)
kw_plot_spec_2D = update_dict(plot_spec_2D_kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)

# get the rvals we want
radlevs = get_slice_levels(dirname)
irvals = kw.irvals
if not kw.rvals is None: # irvals haven't been set directly
    if np.all(kw.rvals == 'all'):
        irvals = np.arange(radlevs.nr)
    else:
        irvals = np.zeros_like(kw.rvals, dtype='int')
        for i in range(len(kw.rvals)):
            irvals[i] = np.argmin(np.abs(radlevs.radius/rsun - kw.rvals[i]))

# and the qvals
qvals = make_array(kw.qvals)

# everything must be array
irvals = make_array(irvals)
rvals = radlevs.radius[irvals]/rsun

modes = make_array(kw.modes)

print (buff_line)
if kw.the_file is None:
    print ("plotting temporal spectra")
    print ("qvals = ", qvals)
    print ("irvals = ", irvals)
    nfigures = len(irvals)*len(qvals)*len(modes)
else:
    print ("plotting ", kw.the_file)
    nfigures = len(modes)
print ("modes = ", modes)
print ("nfigures = ", nfigures)
print (buff_line)

# get and make plotdir of non existent
plotdir = my_mkdir(clas0['plotdir'] + 'tspec/')

for qval in qvals:
    for irval in irvals:
        # get radial level
        rval = rvals[irval]

        # get data
        if kw.the_file is None:
            dataname = ('tspec_qval%04i_irval%02i' %(qval, irval)) + clas0['tag']
           
            kw.the_file = get_widest_range_file(clas0['datadir'], dataname)
        else:
            dataname = get_dataname_from_file(kw.the_file)
        iter1, iter2 = get_iters_from_file(kw.the_file)

        print ("plotting: " + kw.the_file)
        di = get_dict(kw.the_file)
        freq = di['freq']
        vals = np.real(di['vals'])
        nfreq, nell, nm = np.shape(vals)

        for mode in modes:
            basename = mode
            if mode == 'lpower':
                xlabel = 'sph. harm. deg. (l)'
                x = np.arange(nell)
                power = np.sum(vals, axis=2)
                power = power.T/nm # keep power normalized
            if mode == 'mpower':
                xlabel = 'azimuthal wavenumber (m)'
                x = np.arange(nm)
                power = np.sum(vals, axis=1)
                power = power.T/nell # keep power normalized
            if mode == 'lval':
                xlabel = 'azimuthal wavenumber (m)'
                x = np.arange(nm)
                lval = int(kw.lval)
                basename += "%03i" %lval
                power = vals[:, lval, :].T
            if mode == 'mval':
                xlabel = 'sph. harm. deg. (l)'
                x = np.arange(nell)
                mval = int(kw.mval)
                basename += "%03i" %mval
                power = vals[:, :, mval].T
            if mode == 'sect':
                xlabel = 'sectoral (l = m) degree'
                x = np.arange(nell)
                power = np.zeros((nell, nfreq))
                for il in range(nell):
                    power[il, :] = vals[:, il, il]

            # Display at terminal what we are plotting
             
            savename = dataname + '_' + basename + '_' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8)

            # make plot
            fig, axs, fpar = make_figure(**kw_make_figure)
            ax = axs[0, 0]

            #kw_plot_spec_2D.cbar_pos = 'right'
            if kw_plot_spec_2D.x is None:
                kw_plot_spec_2D.x = x
            if kw_plot_spec_2D.y is None:
                kw_plot_spec_2D.y = freq

            plot_spec_2D(power, fig, ax, **kw_plot_spec_2D)

            # make labels
            ylabel = 'freq (Hz)'
            ax.set_xlabel(xlabel, fontsize=default_labelsize)
            ax.set_ylabel(ylabel, fontsize=default_labelsize)

            # make title
            time_string = get_time_string(dirname, iter1, iter2, oneline=True)
            slice_info = ('qval = %04i' %qval) + 5*' ' + (r'$r/R_\odot\ =\ %0.3f$' %rval)

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
