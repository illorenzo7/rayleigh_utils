import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
sys.path.append(os.environ['rapl'] + '/slice')
from common import *
from plotcommon import *
from cla_util import *
from slice_util import spec_2D_fig_dimensions

# Get CLAs
args = sys.argv 
clas0, clas = read_clas(args)
dirname = clas0.dirname
dirname_stripped = strip_dirname(dirname)

# Get the Rayleigh data directory
radatadir = dirname + '/Shell_Slices/'

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists_all(radatadir)

# read first file for some metadata
a0 = Shell_Slices(radatadir + file_list[0], '')

# SPECIFIC ARGS
kwargs_default = dotdict(dict({'the_file': None, 'irvals': np.array([0]), 'rvals': None, 'qvals': None, 'modes': [], 'latvals': [], 'mvals': [], 'diffrot': False, 'azfile': None, 'mmax': None, 'mnot0': False, 'imag': False, 'abs': False}))
# "modes" can be: latpower, mpower, or combinations thereof
# can also get other modes via --latvals (plot freq vs m), --mvals (freq vs lat.)

# many kwargs
kwargs_default.update(my_pcolormesh_kwargs_default)

spec_2D_fig_dimensions['sub_margin_top_inches'] = 1
make_figure_kwargs_default.update(spec_2D_fig_dimensions)
kwargs_default.update(make_figure_kwargs_default)

find_bad_keys(kwargs_default, clas, 'plot/timetrace/tmspec', justwarn=True)

kw = update_dict(kwargs_default, clas)
kw_my_pcolormesh = update_dict(my_pcolormesh_kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)

# check if we want the real or imaginary vals
if kw.imag:
    part = 'imag'
elif kw.abs:
    part = 'abs'
    kw_my_pcolormesh.posdef = True
    kw_my_pcolormesh.logscale = True
else:
    part = 'real'

# get the rvals we want
radlevs = get_slice_levels(dirname)
irvals = kw.irvals

# mmax (if needed)
mmax = kw.mmax

# get latitude
gi = get_grid_info(dirname)
tt_lat = gi['tt_lat']
tw = gi['tw']
tw_nd = tw.reshape((1, 1, len(tw)))

# everything must be an array
irvals = make_array(irvals)
kw.rvals = make_array(kw.rvals)
if not kw.rvals is None: # irvals haven't been set directly
    if isall(kw.rvals):
        irvals = np.arange(radlevs.nr)
    else:
        irvals = np.zeros_like(kw.rvals, dtype='int')
        for i in range(len(kw.rvals)):
            irvals[i] = np.argmin(np.abs(radlevs.radius/rsun - kw.rvals[i]))

# and the qvals
qvals = kw.qvals # ... if groupname is specified, this will just be
                 # the qvals associated with the group, e.g., 
                 # b <--> 801, 802, 803

if isall(qvals): # probably won't use this option here ... would 
    # make too many panels
    qvals = np.sort(a0.qv)

if qvals is None: 
    qvals = np.array([1])

# modes, etc.; these must be lists
modes = make_array(kw.modes, tolist=True)
latvals = make_array(kw.latvals, tolist=True)
mvals = make_array(kw.mvals, tolist=True)

everything = modes + latvals + mvals
if len(everything) == 0:
    modes = ['mpower']
    everything = modes

print (buff_line)
if kw.the_file is None:
    print ("plotting temporal/m spectra")
    print ("qvals = ", qvals)
    print ("irvals = ", irvals)
    nfigures = len(irvals)*len(qvals)*(len(modes) + len(latvals) + len(mvals))
else:
    print ("plotting ", kw.the_file)
    nfigures = len(modes) + len(latvals) + len(mvals)
print ("modes = ", modes)
print ("latvals = ", latvals)
print ("mvals = ", mvals)
print ("nfigures = ", nfigures)
print("considering part %s of the complex spectra" %part)
print (buff_line)

# get and make plotdir if non-existent
if mmax is None:
    plotdir = my_mkdir(clas0['plotdir'] + 'tmspec' + clas0['tag'] + '/')
else:
    plotdir = my_mkdir(clas0['plotdir'] +\
        'tmspec_mmax%03i' %mmax + clas0['tag'] + '/')

# include differential rotation if requested
if kw.diffrot:
    if kw.azfile is None:
        kw.azfile = get_widest_range_file(clas0['datadir'], 'AZ_Avgs')
    print ('over plotting differential rotation')
    print ('reading ' + kw.azfile)
    print (buff_line)
    azdi = get_dict(kw.azfile)
    azvals = azdi['vals']
    azlut = azdi['lut']
    xx = gi['xx']
    om = azvals[:, :, azlut[3]]/xx
    om /= 2*np.pi # rad/s --> Hz
 
for qval in qvals:
    for irval in irvals:
        # get radial level
        rval = radlevs.radius[irval]/rsun
        if kw.diffrot:
            iirval = radlevs.inds[irval]
            # lat profile of diffrot
            om_lat = om[:, iirval]

        # get data
        if kw.the_file is None:
            dataname = ('tmspec_qval%04i_irval%02i' %(qval, irval)) 
            if mmax is None:
                datadir = clas0['datadir'] + 'tmspec/'
            else:
                datadir = clas0['datadir'] + 'tmspec_mmax%03i/' %mmax
            the_file = get_widest_range_file(datadir, dataname)
        else:
            dataname = get_dataname_from_file(kw.the_file)
            the_file = kw.the_file
        iter1, iter2 = get_iters_from_file(the_file)

        print ("reading " + the_file)
        di = get_dict(the_file)
        freq = di['freq']

        vals = di['vals']
        if part == 'imag':
            vals = np.imag(vals)
        elif part == 'abs':
            vals = np.abs(vals)**2
        else:
            vals = np.real(vals)

        # everything with m >= 1 should be counted twice,
        # if looking at power
        if part == 'abs':
            vals[:, 1:, :] *= 2.

            # may want to ignore m = 0 
            # (again, only relevant if we are summing over m, in 
            # modes like "latpower"
            if kw.mnot0:
                vals[:, 0, :] = 0.0

        nfreq, nm, nt = np.shape(vals)

        # add mvals/latvals to modes

        count = 0
        for mode in everything:
            if count < len(modes):
                basename = mode
                if mode == 'latpower':
                    xlabel = 'latitude (deg)'
                    x = tt_lat
                    power = np.sum(vals, axis=1)
                    power = power.T/nm # keep power normalized
                    if kw.diffrot:
                        omy = om_lat
                if mode == 'mpower':
                    xlabel = 'azimuthal wavenumber (m)'
                    x = di['mvals']
                    power = np.sum(vals*tw_nd, axis=2)
                    power = power.T/len(tt_lat) # keep power normalized
                    if kw.diffrot:
                        #print (om_lat)
                        #print (tw)
                        omy = np.sum(om_lat*tw)
                        #print(omy)
                        omy = omy*x
            elif count < len(modes) + len(latvals):
                latval = float(mode)
                ithval = np.argmin(np.abs(tt_lat - latval))
                basename = "ithval%03i" %ithval
                xlabel = 'azimuthal wavenumber (m)'
                x = di['mvals']
                power = vals[:, :, ithval].T
                if kw.diffrot:
                    omy = om_lat[ithval]
                    omy = omy*x
            elif count < len(modes) + len(latvals) + len(mvals):
                mval = int(mode)
                basename = "mval%03i" %mval
                xlabel = 'latitude (deg)'
                x = tt_lat
                power = vals[:, mval, :].T
                if kw.diffrot:
                    omy = om_lat*mval

            if kw.mnot0:
                basename += '_mnot0'

            # Display at terminal what we are plotting
            savename = dataname + '_' + basename + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '_' + part + '.png'

            # make plot
            fig, axs, fpar = make_figure(**kw_make_figure)
            ax = axs[0, 0]

            kw_my_pcolormesh.cbar_pos = 'right'
            if kw_my_pcolormesh.x is None:
                kw_my_pcolormesh.x = x
            if kw_my_pcolormesh.y is None:
                kw_my_pcolormesh.y = freq

            themin, themax = my_pcolormesh(power, fig, ax, **kw_my_pcolormesh)

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
            title += '\nminmax = %1.3e, %1.3e' %(themin, themax)
            title += '\n%s part' %part
            ax.set_title(title, va='bottom', fontsize=default_titlesize)

            # possibly overplot DR
            if kw.diffrot:
                # positive diffrot means negative freq
                omy *= -1
                #print (omy)
                #print (om)
                #print (x)
                # make sure to reset xlim ylim after
                xmin, xmax = ax.get_xlim()
                ymin, ymax = ax.get_ylim()
                ax.scatter(x, omy, color='k')
                ax.set_xlim(xmin,xmax)
                ax.set_ylim(ymin,ymax)


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
