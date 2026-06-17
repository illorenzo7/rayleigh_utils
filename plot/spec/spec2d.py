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

# set fontsize
fontsize = 14

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0.dirname
dirname_stripped = strip_dirname(dirname)

# SPECIFIC ARGS
kw_default = dotdict(dict({'irvals': 0, 'rvals': None, 'rav': True, 'varname': 'vr', 'groupname': None, 'prepend': False, 'dpi': 300, 'minmaxlm': None, 'minmaxl': None, 'minmaxm': None, 'minmax': None, 'log': True}))

# fig dimensions
# 1 row of 3 figures side by side: Power vs lm, Power vs l, Power vs m
fig_dimensions = dotdict({'width_inches': 10., 'sub_aspect': 1, 'sub_margin_right_inches': 1/2, 'nrow': 1, 'ncol': 3})

# Rayleigh data dir
radatadir = dirname + '/Shell_Spectra/'

# now we can update the default kw
kw_my_pcolormesh_default['cbar_pos'] = 'right'
kw_my_pcolormesh_default['posdef'] = True
kw_default.update(kw_my_pcolormesh_default)
kw_make_figure_default.update(fig_dimensions)
kw_default.update(kw_make_figure_default)
for key in range_options: # add the range options key
    kw_default[key] = None
find_bad_keys(kw_default, clas, 'plot/spec/spec2d.py', justwarn=True)

# update relevant keyword args
kw = update_dict(kw_default, clas)
kw_make_figure = update_dict(kw_make_figure_default, clas)
kw_my_pcolormesh = update_dict(kw_my_pcolormesh_default, clas)

# get the data
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], 'Shell_Spectra')
print ('Getting data from ' + kw.the_file)
di = get_dict(kw.the_file)
vals = di['vals']
lut = di['lut']

# get desired sampling locations
sliceinfo = get_sliceinfo(dirname, 'Shell_Spectra')
if kw.rav: # just average over radial direction (simple)
    vals = np.mean(vals, axis=2, keepdims=True)
    nrvals = 1
else: # plotting specific level
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

if kw.groupname is None: # it's a simple variable
    varlabel = kw.varname
    print ("plotting varname = " + kw.varname)
    qval = var_indices[kw.varname]
    field_all_radii = vals[..., lut[qval]]
else:
    varlabel = kw.groupname
    print ("plotting sum(" + kw.varname + ")^2")
    the_group = get_quantity_group(kw.groupname, clas0['magnetism'])
    qvals = the_group['qvals']
    field_all_radii = np.zeros_like(vals[...,0])
    for qval in qvals:
        field_all_radii += vals[..., lut[qval]]

# make the plotting directory
iter1, iter2 = get_iters_from_file(kw.the_file)
plotdir = clas0['plotdir'] + '/spec2d_' + clas0['tag']
plotdir = my_mkdir(plotdir)

# prepare the grid
nell = nm = vals.shape[0]
lvals = np.arange(nell)
mvals = np.arange(nm)
xx, yy = np.meshgrid(lvals, mvals, indexing='ij')

# loop over r-values and make plots
for irval in range(nrvals):
    rval = sliceinfo.rvals[irval]
    savename = 'spec2d_' + str(iter1).zfill(8) +\
            '_' + str(iter2).zfill(8) + '_' + varlabel

    # figure out the r-label
    if kw.rav:
        rlabel = 'rav'
    else:
        rlabel = 'rval%.3f' %rval

    savename += '_' + rlabel + '.png'

    if kw.prepend:
        savename = dirname_stripped + '_' + savename

    savefile = plotdir + '/' + savename
    print("plotting...")

    # now make the plot
    fig, axs, fpar = make_figure(**kw_make_figure)
    field = field_all_radii[..., irval]

    # make 2D power spectrum
    ax = axs[0, 0]
    my_pcolormesh(field, fig, ax, **kw_my_pcolormesh)
    ax.set_xlabel(r'$\ell$', fontsize=fontsize)
    ax.set_ylabel(r'$m$', fontsize=fontsize)
    ax.set_title('2D power spectrum', fontsize=fontsize)

    # make power spectrum vs ell
    ax = axs[0, 1]
    spec_l = np.sum(field, axis=1)
    ax.bar(lvals, spec_l, width=1)

    print ("saving", savefile)
    plt.savefig(savefile, dpi=300)

# close figure at end of loop
if nrvals == 1:
    plt.show()
plt.close()
