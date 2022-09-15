##################################################################
# Routine to plot arbitrary variables (AZ_Avgs)
# Author: Loren Matilsky
# Created: 01/28/2019
##################################################################
# This script plots the quantitities specified by --qvals
# default is v
##################################################################

import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
from azav_util import *
from common import *
from plotcommon import *
from derived_quantities import *
from cla_util import *

# Read command-line arguments (CLAs)
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)
# See if magnetism is "on"
magnetism = clas0['magnetism']

# defaults
kwargs_default = dict({'the_file': None, 'qvals': None, 'groupname': 'v', 'shearfile': None, 'subcurv': False})
kwargs_default.update(plot_azav_grid_kwargs_default)

# overwrite defaults
kw = update_dict(kwargs_default, clas)
kw_plot_azav_grid = update_dict(plot_azav_grid_kwargs_default, clas)

# deal with desired quantities
if kw.qvals is None: # it's a quantity group
    qgroup = get_quantity_group(kw.groupname, magnetism)
    kw.qvals = qgroup['qvals']
    kw_plot_azav_grid.titles = qgroup['titles']
    kw_plot_azav_grid.totsig = qgroup['totsig']
    kw_plot_azav_grid.ncol = qgroup['ncol']
else:
    kw_plot_azav_grid.titles = parse_quantities(kw.qvals)[1]
    kw.groupname = input("choose a groupname to save your plot\n to not save it, enter 'nosave': ")
# check for bad keys
find_bad_keys(kwargs_default, clas, clas0['routinename'], justwarn=True)

# get the data type we want; not all will be an AZ_Avgs file
dataname_list = dict({})
for ext in ['tot', 'pmp', 'ppm', 'mmm', 'mpp', 'ppp']:
    for ext2 in ['', 'r', 't', 'p']:
        dataname_list['meprodnum' + ext + ext2] = 'me_prod'
        dataname_list['meprodshear' + ext + ext2] = 'me_prod_shear'
        dataname_list['meprodadvec' + ext + ext2] = 'me_prod_advec'
for direc in ['r', 't', 'p']:
    dataname_list['ind' + direc + 'alt'] = 'induct_alt'
    dataname_list['ind' + direc + 'altnum'] = 'induct_alt_num'
for ext in ['mm', 'ms']:
    dataname_list['magtorque' + ext] = 'mag_torque'
for ext in ['tot', 'mmm', 'mpp']:
    dataname_list['meprodmean' + ext] = 'me_prod_mean'
dataname_list['ferraro'] = 'ferraro'
dataname_list['ferraroinst'] = 'ferraroinst'

if kw.groupname in dataname_list.keys():
    dataname = dataname_list[kw.groupname]
else:
    dataname = 'AZ_Avgs'

# get data
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], dataname)

print ('Getting data from ' + kw.the_file)
di = get_dict(kw.the_file)
vals = di['vals']
if dataname == 'AZ_Avgs':
    lut = di['lut']

if kw.subcurv: # want to subtract curvature terms from exact ME prod terms
    # only makes sense for exact me prod, or me_prod_[shear, advec]
    print ('subcurv = True')
    if 'meprod' in [kw.groupname[:-3], kw.groupname[:-4]]: # this is the exact stuff
        # should only specify subcurv for exact quantity groups
        baselen = 6
    elif kw.groupname[:11] in ['meprodshear', 'meprodadvec']:
        baselen = 11
    else:
        print ('ERROR: subcurv = True meaningless for groupname = ' + kw.groupname)
        print ('exiting')
        sys.exit()

    # get data type (meprod, meprodshear, meprodadvec)
    datatype = kw.groupname[:baselen]
    ext = kw.groupname[baselen:baselen + 3] # tot, pmp, etc.

    if datatype == 'meprod': # need curvature terms directly
        print ('subtracting canceling curvature terms from ' + kw.groupname)
        if kw.shearfile is None:
            kw.shearfile = get_widest_range_file(clas0['datadir'], 'me_prod_shear')
        print ('reading curvature terms from ' + kw.shearfile)
        di_shear = get_dict(kw.shearfile)
        vals_shear = di_shear['vals']
        if ext == 'tot':
            iqshear = 0
        if ext == 'pmp':
            iqshear = 15
        if ext == 'ppm':
            iqshear = 30
        if ext == 'mmm':
            iqshear = 45
        if ext == 'mpp':
            iqshear = 60
        if ext == 'ppp':
            iqshear = 75
        curvr1 = vals_shear[..., iqshear + 3]
        curvr2 = vals_shear[..., iqshear + 4]
        curvt2 = vals_shear[..., iqshear + 9]

# collect terms to plot
print ("plotting the following quantities:")
print ("qvals = ", kw.qvals)

terms = []
for qval in kw.qvals:
    if is_an_int(qval):
        qval = int(qval)
        if dataname == 'AZ_Avgs':
            if lut[qval] == 4000: # derived quantity
                terms.append(derive_quantity(dirname, vals, lut, qval)[..., 0])
            else:
                terms.append(vals[:, :, lut[qval]])
        else:
            terms.append(vals[:, :, qval])
    else:
        terms.append(derived_azav(dirname, vals, lut, qval))

if kw.subcurv:
    # first need to figure out if this is only for a specific direction
    ext2 = kw.groupname[baselen+3:]

    if datatype == 'meprod':
        if ext2 == 'r':
            # radial quantities
            terms[2] -= (curvr1 + curvr2) # subtract from shear
            terms[3] += (curvr1 + curvr2) # add to advection
        elif ext2 == 't':
            terms[2] -= curvt2
            terms[3] += curvt2
        elif ext2 == '':
            # radial quantities
            terms[2] -= (curvr1 + curvr2) # subtract from shear
            terms[3] += (curvr1 + curvr2) # add to advection
            # theta quantities
            terms[2 + qgroup['ncol']] -= curvt2
            terms[3 + qgroup['ncol']] += curvt2
        else:
            # nothing to do for phi quantities
            print ('WARNING: subcurv = True has no effect for qgroup ' + kw.groupname)
    else: # meprodshear, meprodadvec
        print ('got here')
        if ext2 == 'r':
            kw_plot_azav_grid.totsig[3:5] = 0
        elif ext2 == 't':
            kw_plot_azav_grid.totsig[4] = 0
        elif ext2 == '':
            kw_plot_azav_grid.totsig = np.ones(len(terms))
            kw_plot_azav_grid.totsig[3:5] = 0
            kw_plot_azav_grid.totsig[9] = 0
        else:
            # nothing to do for phi quantities
            print ('WARNING: subcurv = True has no effect for qgroup ' + kw.groupname)

# make the main title
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2)

if kw.groupname is None:
    qlabel = array_of_strings(kw.qvals)
else:
    qlabel = kw.groupname

if kw_plot_azav_grid.maintitle is None:
    kw_plot_azav_grid.maintitle = dirname_stripped + '\n' + 'qvals = ' + qlabel + '\n' + time_string

# Generate the figure using standard routine
di_grid = get_grid_info(dirname)
if kw.shav:
    kw_plot_azav_grid.tw = di_grid['tw']
figs = plot_azav_grid (terms, di_grid['rr'], di_grid['cost'], **kw_plot_azav_grid)

if kw.shav:
    fig, av_fig = figs
else:
    fig = figs

# save the figure if tag (or qgroup) was specified
if len(clas0['tag']) > 0 or not kw.groupname is None:
    basename = 'azav_'
    if not kw.groupname is None:
        basename += kw.groupname
    basename += clas0['tag']

if basename in ['azav_v', 'azav_b']: # these go in main directory
    plotdir = my_mkdir(clas0['plotdir'])
else:
    plotdir = my_mkdir(clas0['plotdir'] + 'azav/')

if clas0['saveplot'] and kw.groupname != 'nosave':
    savefile = plotdir + basename + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
    print ('saving figure at ' + savefile)
    fig.savefig(savefile, dpi=300)

    if kw.shav:
        basename = basename.replace('azav', 'shav')
        av_savefile = plotdir + basename + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
        print ('saving lat. avg. figure at ' + av_savefile)
        av_fig.savefig(av_savefile, dpi=300)

if clas0['showplot']:
    plt.show()

plt.close(fig)
if kw.shav:
    plt.close(av_fig)
