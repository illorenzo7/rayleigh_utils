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
kwargs_default = dotdict(dict({'the_file': None, 'irvals': np.array([0]), 'rvals': None, 'qvals': np.array([1]), 'modes': [], 'lvals': [], 'mvals': [], 'sectoff': []}))
# "modes" can be: lpower, mpower, sect, or combinations thereof
# can also get other modes via --lvals (plot freq vs m), --mvals (freq vs l), and --sectoff (m = l - offset)

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

# everything must be an array
irvals = make_array(irvals)

# modes, etc.; these must be lists
modes = make_array(kw.modes, tolist=True)
lvals = make_array(kw.lvals, tolist=True)
mvals = make_array(kw.mvals, tolist=True)
sectoff = make_array(kw.sectoff, tolist=True)

everything = modes + lvals + mvals + sectoff
if len(everything) == 0:
    modes = ['lpower']
    everything = modes

print (buff_line)
if kw.the_file is None:
    print ("plotting temporal spectra")
    print ("qvals = ", qvals)
    print ("irvals = ", irvals)
    nfigures = len(irvals)*len(qvals)*(len(modes) + len(lvals) + len(mvals) + len(sectoff))
else:
    print ("plotting ", kw.the_file)
    nfigures = len(modes) + len(lvals) + len(mvals) + len(sectoff)
print ("modes = ", modes)
print ("lvals = ", lvals)
print ("mvals = ", mvals)
print ("sectoff = ", sectoff)
print ("nfigures = ", nfigures)
print (buff_line)

# get and make plotdir of non existent
plotdir = my_mkdir(clas0['plotdir'] + 'tspec/')
 
for qval in qvals:
    for irval in irvals:
        # get radial level
        rval = radlevs.radius[irval]/rsun

        # get data
        if kw.the_file is None:
            dataname = ('tspec_qval%04i_irval%02i' %(qval, irval))          
            the_file = get_widest_range_file(clas0['datadir'], dataname)
        else:
            dataname = get_dataname_from_file(kw.the_file)
            the_file = kw.the_file
        iter1, iter2 = get_iters_from_file(the_file)

        print ("reading " + the_file)
        di = get_dict(the_file)
        freq = di['freq']
        vals = np.real(di['vals'])
        nfreq, nell, nm = np.shape(vals)

        # add mvals/lvals sectoff to modes

        count = 0
        for mode in everything:
            if count < len(modes):
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
                if mode == 'sect':
                    xlabel = 'sectoral (l = m) degree'
                    x = np.arange(nell)
                    power = np.zeros((nell, nfreq))
                    for il in range(nell):
                        power[il, :] = vals[:, il, il]
            elif count < len(modes) + len(lvals):
                lval = int(mode)
                basename = "lval%03i" %lval
                xlabel = 'azimuthal wavenumber (m)'
                x = np.arange(nm)
                power = vals[:, lval, :].T
            elif count < len(modes) + len(lvals) + len(mvals):
                mval = int(mode)
                basename = "mval%03i" %mval
                xlabel = 'sph. harm. deg. (l)'
                x = np.arange(nell)
                power = vals[:, :, mval].T
            elif count < len(modes) + len(lvals) + len(mvals) + len(sectoff):
                basename = "sectoff%03i" %mode
                offset = int(mode)
                xlabel = 'sectoral offset (m = l - %i)' %offset
                x = np.arange(nell)
                power = np.zeros((nell, nfreq))
                for il in range(offset, nell):
                    power[il, :] = vals[:, il, il - offset]

            # Display at terminal what we are plotting
            savename = dataname + '_' + basename + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

            # make plot
            fig, axs, fpar = make_figure(**kw_make_figure)
            ax = axs[0, 0]

            kw_plot_spec_2D.cbar_pos = 'right'
            if kw_plot_spec_2D.y is None:
                kw_plot_spec_2D.y = freq

            plot_spec_2D(power, fig, ax, **kw_plot_spec_2D)

            # make labels
            ylabel = 'freq (Hz)'
            ax.set_xlabel(xlabel, fontsize=default_labelsize)
            ax.set_ylabel(ylabel, fontsize=default_labelsize)

            # get y axis to use scientific notation
            ax.ticklabel_format(scilimits=(-3,4), useMathText=True)

            # make title
            time_string = get_time_string(dirname, iter1, iter2, oneline=True)
            slice_info = ('qval = %04i' %qval) + 5*' ' + (r'$r/R_\odot\ =\ %0.3f$' %rval) + 5*' ' + "mode = " + basename

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
            count += 1
print (buff_line)
