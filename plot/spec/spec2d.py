# Author: Loren Matilsky
# Date created: 11/29/2022
# Description: go-to slice plotting routine for everything

# import modules
import sys, os
sys.path.append(os.environ['raco'])
from common import *
import matplotlib.pyplot as plt
import numpy as np
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['rapl'])
from plotcommon import *
from cla_util import *
#from get_slice import get_slice, get_label

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0.dirname
dirname_stripped = strip_dirname(dirname)

# SPECIFIC ARGS
kw_default = dotdict(dict({'irvals': 0, 'rvals': None, 'varnames': 'vr', 'prepend': False, 'dpi': 300, 'minmaxl': None, 'minmaxm': None, 'log': False}))

# fig dimensions
# 1 row of 3 figures side by side: Power vs lm, Power vs l, Power vs m
fig_dimensions = dotdict({'width_inches': 10., 'sub_aspect': 1, 'sub_margin_right_inches': 1/2, 'nrow': 1, 'ncol': 3})
rlabel = 'rval'
rfmt = '%1.3e'

# Rayleigh data dir
radatadir = dirname + '/Shell_Spectra/'

# now we can update the default kw
kw_default.update(kw_plotting_func_default)
kw_default.update(kw_my_pcolormesh_default)
kw_make_figure_default.update(fig_dimensions)
kw_default.update(kw_make_figure_default)
for key in range_options: # add the range options key
    kw_default[key] = None
find_bad_keys(kw_default, clas, 'plot/spec/spec2d.py', justwarn=True)

# update relevant keyword args
kw = update_dict(kw_default, clas)
kw_plotting_func = update_dict(kw_plotting_func_default, clas)
kw_make_figure = update_dict(kw_make_figure_default, clas)
kw_my_pcolormesh = update_dict(kw_my_pcolormesh_default, clas)

# get the data
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], 'Shell_Spectra')
print ('Getting data from ' + kw.the_file)
di = get_dict(kw.the_file)
vals = di['vals']
lut = di['lut']

# get desired varnames
sliceinfo = get_sliceinfo(dirname, 'Shell_Spectra')
if isall(kw.varnames):
    kw.varnames = array_of_strings(sliceinfo.qv)
# no matter what, this needs to be an array of strings
kw.varnames = array_of_strings(make_array(kw.varnames))
nq = len(kw.varnames)

# get desired sampling locations
# first the indices must be an array
kw.irvals = make_array(kw.irvals)

# get the rvals we want 
if not kw.rvals is None: # rvals have been set directly
    # need the available sampling locations
    if isall(kw.rvals):
        kw.irvals = np.arange(sliceinfo.nrvals)
    else:
        kw.rvals = make_array(kw.rvals)
        kw.irvals = inds_from_vals(sliceinfo.rvals, kw.rvals)

# these are the r vals we end up with
nrvals = len(kw.rvals)

# say what we are plotting
print (buff_line)

# print varnames
print (buff_line)
print (("plotting %i variables:\nvarnames = " %nq) +\
        arr_to_str(kw.varnames, "%s"))

# calculate total number of figures
nfigures = nrvals*nq
print (buff_line)
print ("nfigures = %i x %i = %i"\
        %(nrvals, nq, nfigures))
print (buff_line)

for varname in kw.varnames:
    varlabel = get_label(varname)
    simple_label = varname
    for irval in kw.irvals:
        rval = sliceinfo.rvals[irval]
        # now we're in the big loop
        # get plot directory and image name to save the file
        plotdir = clas0['plotdir'] + plottype + clas0['tag']
        savename = 'spec2d_' + fname + '_' + simple_label +\
            (('_' + rlabel + rfmt) %rval) + '.png'
        if kw.prepend:
            savename = dirname_stripped + '_' + savename

        savefile = plotdir + '/' + savename
        plotdir = my_mkdir(plotdir)

    # get the slice
    if plottype in ['moll', 'ortho']:
        varname_root, deriv, primevar, sphvar = get_varprops(varname)
        qvals = None # by default, but update
        if varname_root in var_indices:
            qvals = var_indices[varname_root]
        elif varname_root == 'omz':
            qvals = [301, 302]
        a = reading_func(radatadir + fname, '', irvals=irval, qvals=qvals)
    else:
        a = reading_func(radatadir + fname, '')

    # get the variable
    vals = get_slice(a, varname, dirname=dirname)

    # get r location
    if dataname == 'Shell_Slices':
        field = vals[:, :, 0] # only read in one r val
    elif dataname == 'Meridional_Slices':
        field = vals[irval, :, :]
    else: # equatorial slices
        field = vals

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

    # possibly add lon avg.
    if lonav:
        fig_width_inches, fig_height_inches = fig.get_size_inches()
        fig_aspect = fig_height_inches/fig_width_inches
        # get ax dimensions
        ax_left, ax_right, ax_bottom, ax_top = axis_range(ax)
        ax_width = ax_right - ax_left
        ax_height = ax_top - ax_bottom

        ax_line_left = ax_right + 1/8/fig_width_inches
        ax_line_bottom = ax_bottom
        ax_line_width = 1 - 5/8/fig_width_inches - ax_line_left
        ax_line_height = ax_height
        ax_line = fig.add_axes([ax_line_left, ax_line_bottom, ax_line_width, ax_line_height])

        # make line plot
        field_lonav = np.mean(field, axis=0)
        tt_lat = (np.pi/2. - np.arccos(a.costheta))*180./np.pi
        ax_line.plot(field_lonav, tt_lat)

        # set ylim
        ax_line.set_ylim(-90., 90.)

        # mark zero lines
        mark_axis_vals(ax_line, 'x')
        mark_axis_vals(ax_line, 'y')

        # move ticks to right
        ax_line.yaxis.tick_right()
        ax_line.yaxis.set_label_position("right")

        # Get ticks everywhere
        plt.sca(ax_line)
        plt.minorticks_on()
        plt.tick_params(top=True, right=True, direction='in', which='both')

        # label the axes
        ax_line.set_title('lon. avg.')
        ax_line.set_ylabel('latitude (degrees)')

    time_string = get_time_string(dirname, t1=a.time[0],SF=5)

    if plottype == 'moll':
        location_and_perspective =\
            (rlabel + ' = ' + rfmt) %rval +\
                ('\n' + r'$\phi_0$' + ' = ' + lon_fmt_tex) %clon
    elif plottype == 'ortho':
        location_and_perspective =\
            (rlabel + ' = ' + rfmt) %rval +\
                ('\n' + r'$\phi_0$' + ' = ' + lon_fmt_tex) %clon +\
                ('    ' + r'$\lambda_0$' + ' = ' + lat_fmt_tex) %clat
    elif plottype == 'eq':
        location_and_perspective = (r'$\phi_0$' + ' = ' + lon_fmt_tex) %clon
    elif plottype == 'mer':
        location_and_perspective = (rlabel + ' = ' + rfmt) %rval

    title = dirname_stripped + '\n' +\
            varlabel + '\n' +\
            location_and_perspective + '\n' +\
            time_string
    ax.set_title(title, va='bottom', fontsize=default_titlesize)

    # save by default
    if clas0['saveplot']:
        plt.savefig(savefile, dpi=kw.dpi)
    if rank == 0:
        pcnt_done = (ifigure+1)/my_nfigures*100.
        # print what we saved and how far along we are
        if clas0['saveplot']:
            print(fill_str("saved " + savefile + " (dpi=%i)" %kw.dpi) + '\n' +\
                    ('rank 0 %5.1f%% done' %pcnt_done), end='\r')
        # always show if nfigures is 1
        if nfigures == 1 and clas0['showplot']:
            plt.show()   
    # close figure at end of loop
    plt.close()

# Checkpoint and time
comm.Barrier()
if rank == 0:
    t2 = time.time()
    print('\n' + fill_str('plotting time'), end='')
    print (format_time(t2 - t1))
    print(make_bold(fill_str('total time')), end='')
    print (make_bold(format_time(t2 - t1_glob)))
    print (buff_line)
