import matplotlib.pyplot as plt
import time
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['rapl'])
sys.path.append(os.environ['rapl'] + '/azav')
from common import *
from plotcommon import *
from cla_util import *
from slice_util import *
from rayleigh_diagnostics import Shell_Slices, Equatorial_Slices, Meridional_Slices
from azav_util import plot_azav
#from get_slice import get_slice, get_label

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0.dirname
dirname_stripped = strip_dirname(dirname)

# SPECIFIC ARGS
kwargs_default = dotdict(dict({'type': None, 'isamplevals': np.array([0]), 'samplevals': None, 'varnames': np.array(['vr'])}))

# this guy need to update right away to choose fig dimensions
if clas.type is None:
    plottype = 'moll'
else:
    plottype = clas.type

print (buff_line)
print ("plottype = " + plottype)
if plottype == 'moll':
    fig_dimensions = moll_fig_dimensions
if plottype == 'ortho':
    fig_dimensions = ortho_fig_dimensions
if plottype in ['moll', 'ortho']:
    plotting_func = plot_moll_or_ortho
    plotting_func_kwargs_default = plot_moll_or_ortho_kwargs_default
    if plottype == 'ortho':
        plotting_func_kwargs_default['ortho'] = True
    dataname = 'Shell_Slices'
    reading_func = Shell_Slices
if plottype == 'eq':
    fig_dimensions = eq_fig_dimensions
    plotting_func = plot_eq
    plotting_func_kwargs_default = plot_eq_kwargs_default
    dataname = 'Equatorial_Slices'
    reading_func = Equatorial_Slices
if plottype == 'mer':
    fig_dimensions = mer_fig_dimensions
    plotting_func = plot_mer
    plotting_func_kwargs_default = plot_mer_kwargs_default
    dataname = 'Meridional_Slices'
    reading_func = Meridional_Slices

# now we can update the default kwargs
kwargs_default.update(plotting_func_kwargs_default)
make_figure_kwargs_default.update(fig_dimensions)
kwargs_default.update(make_figure_kwargs_default)
find_bad_keys(kwargs_default, clas, 'plot/slice/raslice', justwarn=True)

# update relevant keyword args
kw = update_dict(kwargs_default, clas)
kw_plotting_func = update_dict(plotting_func_kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)

# these need to be arrays
kw.isamplevals = make_array(kw.isamplevals)
kw.varnames = array_of_strings(make_array(kw.varnames))

# make plot directory if nonexistent
basename = plottype
plotdir = my_mkdir(clas0['plotdir'] + basename + clas0['tag'] + '/')

# Get desired file names in datadir and their integer counterparts
radatadir = dirname + '/' + dataname + '/'
clas_mod = dict({'iter': 'last'})
clas_mod.update(clas)
file_list, int_file_list, nfiles = get_file_lists(radatadir, clas_mod)

# need one of these no matter what
print ("reading " + dataname + '/' + file_list[0])
t1 = time.time()
a0 = reading_func(radatadir + file_list[0], '')
print ("done reading")
print ('plotting %i %s files: %s through %s'\
    %(nfiles, dataname, file_list[0], file_list[-1]))

# need to know which sample values are available
if plottype in ['moll', 'ortho']:
    samplevals_avail = a0.radius
    samplelabel = 'rval'
    samplefmt = '%1.3e'
if plottype == 'mer':
    samplevals_avail = a0.phi*180/np.pi
    samplelabel = 'lon'
    samplefmt = lon_fmt
# for equatorial slices it doesn't matter -- only one equator!

# get the samplevals we want
if not kw.samplevals is None: # samplevals have been set directly
    if kw.samplevals == 'all':
        kw.isamplevals = np.arange(len(samplevals_avail))
    else:
        kw.samplevals = make_array(kw.samplevals)
        kw.isamplevals = np.zeros_like(kw.samplevals, dtype='int')
        for i in range(len(kw.samplevals)):
            kw.isamplevals[i] = np.argmin(np.abs(samplevals_avail - kw.samplevals[i]))

# these are the sample vals we end up with
if not plottype == 'eq':
    kw.samplevals = samplevals_avail[kw.isamplevals]

# get the vars we want
if kw.varnames == 'all':
    kw.varnames = array_of_strings(a0.qv)

# loop over samplevals/vars and make plots
print (buff_line)
if not plottype == 'eq':
    # print sampling locations
    print ("i%ss = " %samplelabel + arr_to_str(kw.isamplevals, '%i'))
    print ("%ss = " %samplelabel + arr_to_str(kw.samplevals, samplefmt))
print ("varnames =", kw.varnames)
print ("nfiles =", nfiles)
nfigures = len(kw.isamplevals)*len(kw.varnames)*nfiles
print ("nfigures =", nfigures)
print (buff_line)

for fname in file_list:
    if fname == file_list[0]:
        a = a0
    else:
        a = reading_func(radatadir + fname, '')
    t2 = time.time()
    print ("took %1.3e sec" %(t2-t1))
    t1 = t2
    for varname in kw.varnames:
        # get the desired field variable
        vals = get_slice(a, varname, dirname=dirname)

        # variable labels
        basic = is_basic(varname)
        if basic:
            varlabel = get_label(varname)
            simple_label = varname
        else:
            varlabel, simple_label = get_label(varname)

        # loop over the sampling locations
        for isampleval in kw.isamplevals:
            if dataname == 'Shell_Slices':
                field = vals[:, :, isampleval]
            elif dataname == 'Meridional_Slices':
                field = vals[isampleval, :, :]
            else: # equatorial slices
                field = vals[:, :]

            # Display at terminal what we are plotting (and saving)
            savename = basename + '_' + str(a.iters[0]).zfill(8) + '_' + simple_label 
            if not plottype == 'eq':
                sampleval = samplevals_avail[isampleval]
                savename += ('_' + samplelabel + samplefmt) %sampleval
                
            savename += '.png'

            # make plot
            fig, axs, fpar = make_figure(**kw_make_figure)
            ax = axs[0, 0]
            if plottype in ['moll', 'ortho']:
                plotting_args = field, a.costheta, fig, ax
            elif plottype == 'mer':
                plotting_args = field, a.radius, a.costheta, fig, ax
            elif plottype == 'eq':
                plotting_args = field, a.radius, fig, ax

            plotting_func(*plotting_args, **kw_plotting_func)

            # make title
            time_string = get_time_string(dirname, a.iters[0])

            if plottype == 'moll':
                location_and_perspective =\
                    (samplelabel + ' = ' + samplefmt) %sampleval +\
                        ('\nclon = ' + lon_fmt) %kw.clon
            elif plottype == 'ortho':
                location_and_perspective =\
                    (samplelabel + ' = ' + samplefmt) %sampleval +\
                        ('\nclon = ' + lon_fmt) %kw.clon +\
                        '   clat = ' + lat_format(kw.clat)
            elif plottype == 'eq':
                location_and_perspective = ('clon = ' + lon_fmt) %kw.clon 
            elif plottype == 'mer':
                location_and_perspective = (samplelabel + ' = ' + samplefmt) %sampleval

            title = dirname_stripped + '\n' +\
                    varlabel + '\n' +\
                    location_and_perspective + '\n' +\
                    time_string
            ax.set_title(title, va='bottom', fontsize=default_titlesize)

            # save by default
            if clas0['saveplot']:
                print ("saving " + plotdir + savename)
                plt.savefig(plotdir + savename, dpi=300)
            # always show if nfigures is 1
            if nfigures == 1 and clas0['showplot']:
                print ("displaying " + plotdir + savename)
                plt.show()   
            t2 = time.time()
            print ("took %1.3e sec" %(t2-t1))
            t1 = t2
            print (buff_line)            
            plt.close()
