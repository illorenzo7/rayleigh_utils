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
fontsize = 12

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0.dirname
dirname_stripped = strip_dirname(dirname)

# SPECIFIC ARGS
kw_default = dotdict(dict({'irvals': 0, 'rvals': None, 'rav': True, 'varname': 'vr', 'groupname': None, 'prepend': False, 'dpi': 300, 'minmaxlm': None, 'minmaxl': None, 'minmaxm': None, 'minmax': None, 'log': False}))

# fig dimensions
# 1 row of 3 figures side by side: Power vs lm, Power vs l, Power vs m
fig_dimensions = dotdict({'width_inches': 10., 'sub_aspect': 1, 'sub_margin_right_inches': 1/2, 'nrow': 2, 'ncol': 3, 'margin_top_inches': 1.25})

# Rayleigh data dir
radatadir = dirname + '/Shell_Spectra/'

# now we can update the default kw
kw_my_pcolormesh_default['cbar_pos'] = 'right'
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
    rvals_nd = sliceinfo.rvals.reshape((1, 1, sliceinfo.nr, 1))
    vals = np.mean(vals*rvals_nd**2, axis=2, keepdims=True)/np.mean(sliceinfo.rvals**2)
    nrvals = 1
else: # plotting specific level
    # first the indices must be an array
    kw.irvals = make_array(kw.irvals)

    # get the rvals we want 
    if not kw.rvals is None: # rvals have been set directly
        # need the available sampling locations
        if isall(kw.rvals):
            kw.irvals = np.arange(sliceinfo.nr)
        else:
            kw.rvals = make_array(kw.rvals)
            kw.irvals = inds_from_vals(sliceinfo.rvals, kw.rvals)

    # these are the r vals we end up with
    nrvals = len(kw.rvals)

# say what we are plotting
print (buff_line)

# print varnames
print (buff_line)
texlabels = dotdict({'v': r'$\mathbf{v}$', 'b': r'$\mathbf{B}$', 'om': r'$\mathbf{\omega}$', 'j': r'$\mathbf{\mathcal{J}}$'})

if kw.groupname is None: # it's a simple variable
    if is_basic(kw.varname):
        varlabel_title = get_label(kw.varname)
        varlabel = kw.varname
    else:
        varlabel_title, varlabel = get_label(kw.varname)
    print ("plotting var = (" + kw.varname + ")^2")
    qval = var_indices[kw.varname]
    field_all_radii = vals[..., lut[qval]]
else:
    varlabel = kw.groupname
    varlabel_title = texlabels[kw.groupname]    
    print ("plotting var = sum(" + kw.groupname + ")^2")
    the_group = get_quantity_group(kw.groupname, clas0['magnetism'])
    qvals = the_group['qvals']
    field_all_radii = np.zeros_like(vals[...,0])
    for qval in qvals:
        field_all_radii += vals[..., lut[qval]]

# make the plotting directory
iter1, iter2 = get_iters_from_file(kw.the_file)
plotdir = clas0['plotdir'] + '/spec2d' + clas0['tag']
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
        rlabel_title = 'volume average'
    else:
        rlabel = 'rval%.3f' %rval
        rlabel_title = '%.3f' %rval

    savename += '_' + rlabel + '.png'

    if kw.prepend:
        savename = dirname_stripped + '_' + savename

    savefile = plotdir + '/' + savename

    # now make the plot
    fig, axs, fpar = make_figure(**kw_make_figure)
    field = field_all_radii[..., irval]

    # make 2D power spectrum
    ax = axs[0, 0]
    my_pcolormesh(field, fig, ax, **kw_my_pcolormesh)
    ax.set_xlabel(r'$\ell$', fontsize=fontsize)
    ax.set_ylabel(r'$m$', fontsize=fontsize)
    the_title = '2D power spectrum of |' + varlabel_title + r'$|^2$'
    ax.set_title(the_title, fontsize=fontsize)

    lmin, lmax = ax.get_xlim()
    mmin, mmax = ax.get_ylim()

    # make power spectrum vs ell
    ax = axs[0, 1]
    spec_l = np.sum(field, axis=1)
    ax.bar(lvals, spec_l, width=1)
    ax.set_xlim(lmin, lmax)
    # set y limits
    if kw.log:
        ymin, ymax = lineplot_minmax(lvals, [spec_l], log=kw.log, buff_ignore=0.05)
        ax.set_ylim(ymin, ymax)
    
    ax.set_xlabel(r'$\ell$', fontsize=fontsize)
    ax.set_ylabel('power', fontsize=fontsize)
    the_title = r'$\ell$' + '-power: ' + r'$\Sigma_{m=0}^\ell|$' + varlabel_title + r'$|^2$'
    ax.set_title(the_title, fontsize=fontsize)

    # make power spectrum vs m
    ax = axs[0, 2]
    spec_m = np.sum(field, axis=0)
    ax.bar(mvals, spec_m, width=1)
    ax.set_xlim(mmin, mmax)
    # set y limits
    if kw.log:
        ymin, ymax = lineplot_minmax(lvals, [spec_m], log=kw.log, buff_ignore=0.05)
        ax.set_ylim(ymin, ymax)

    ax.set_xlabel(r'$m$', fontsize=fontsize)
    ax.set_ylabel('power', fontsize=fontsize)
    the_title = r'$m$' + '-power: ' + r'$\Sigma_{\ell=m}^{\ell_{\rm max}}|$' + varlabel_title + r'$|^2$'

    ax.set_title(the_title, fontsize=fontsize)

    # find ell = m power
    ax = axs[1, 0]
    spec_sectoral = np.diag(field)
    ax.bar(lvals, spec_sectoral, width=1)
    ax.set_xlim(lmin, lmax)
    # set y limits
    if kw.log:
        ymin, ymax = lineplot_minmax(lvals, [spec_sectoral], log=kw.log, buff_ignore=0.05)
        ax.set_ylim(ymin, ymax)
    ax.set_xlabel(r'$\ell=m$', fontsize=fontsize)
    ax.set_ylabel('power', fontsize=fontsize)
    the_title = 'sectoral power: ' + r'$|$' + varlabel_title + r'$|^2(\ell=m)$'
    ax.set_title(the_title, fontsize=fontsize)

    # make sectoral power spectrum vs ell
    field_sectoral_removed = np.copy(field)
    for lmval in lvals:
        field_sectoral_removed[lmval, lmval] = 0.

    ax = axs[1, 1]
    spec_l = np.sum(field_sectoral_removed, axis=1)
    ax.bar(lvals, spec_l, width=1)
    ax.set_xlim(lmin, lmax)
    # set y limits
    if kw.log:
        ymin, ymax = lineplot_minmax(lvals, [spec_l], log=kw.log, buff_ignore=0.05)
        ax.set_ylim(ymin, ymax)
    ax.set_xlabel(r'$\ell$', fontsize=fontsize)
    ax.set_ylabel('power', fontsize=fontsize)
    the_title = 'nonsectoral ' + r'$\ell$' + '-power: ' + r'$\Sigma_{m=0}^{\ell-1}|$' + varlabel_title + r'$|^2$'
    ax.set_title(the_title, fontsize=fontsize)

    # make power spectrum vs m
    ax = axs[1, 2]
    spec_m = np.sum(field_sectoral_removed, axis=0)
    ax.bar(mvals, spec_m, width=1)
    ax.set_xlim(mmin, mmax)
    # set y limits
    if kw.log:
        ymin, ymax = lineplot_minmax(lvals, [spec_m], log=kw.log, buff_ignore=0.05)
        ax.set_ylim(ymin, ymax)
    ax.set_xlabel(r'$m$', fontsize=fontsize)
    ax.set_ylabel('power', fontsize=fontsize)
    the_title = 'nonsectoral ' + r'$m$' + '-power: ' + r'$\Sigma_{\ell=m+1}^{\ell_{\rm max}}|$' + varlabel_title + r'$|^2$'
    ax.set_title(the_title, fontsize=fontsize)

    # possibly set y scale to log
    for ax in axs.flatten()[1:]:
        if kw.log:
            ax.set_yscale('log')

    # now set main title
    t1, t2 = get_times_from_file(dirname, kw.the_file)
    time_string = get_time_string(dirname, t1=t1,t2=t2,SF=5)

    maintitle = dirname_stripped + '\n' +\
            'var = ' + varlabel_title + '\n' +\
            r'$r=$' + rlabel_title + '\n' +\
            time_string
    fig.text(0.5, 1 - 1/8/fpar['height_inches'], maintitle, ha='center', va='top', fontsize=fontsize)

    print ("saving", savefile)
    plt.savefig(savefile, dpi=300)

# close figure at end of loop
if nrvals == 1:
    plt.show()
plt.close()
