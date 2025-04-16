# Author: Loren Matilsky
# Created: 12/19/2022
#
# Description: Script to plot arbitrary list of quantitities 
# (specified by --qvals; default v)

import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
from azav_util import *
from common import *
from plotcommon import *
from cla_util import *

# Read command-line arguments (CLAs)
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)
# See if magnetism is "on"
magnetism = clas0['magnetism']
advect_reference_state = clas0['advect_reference_state']

# defaults
kw_default = dict({'the_file': None, 'qvals': None, 'groupname': None})
kw_default.update(kw_plot_azav_grid_default)

# overwrite defaults
kw = update_dict(kw_default, clas)
kw_plot_azav_grid = update_dict(kw_plot_azav_grid_default, clas)

# deal with possibly different aspect ratios
if kw.halfplane:
    kw_plot_azav_grid.sub_aspect = 1
if kw.modrms:
    kw_plot_azav_grid.sub_aspect += 1
    kw_plot_azav_grid.sub_margin_right_inches += 1/2

# need a bit extra room for subplot labels
kw_plot_azav_grid.sub_margin_top_inches += 1/4

# deal with desired quantities
if kw.qvals is None: # it's a quantity group
    if kw.groupname is None: # by default plot groupname = v
        kw_plot_azav_grid.groupname = kw.groupname = 'v'
    qgroup = get_quantity_group(kw.groupname, magnetism, advect_reference_state)
    kw.qvals = qgroup['qvals']
    kw_plot_azav_grid.titles = list(qgroup['titles'])
    if not kw_plot_azav_grid.totsig is None:
        kw_plot_azav_grid.totsig = list(qgroup['totsig'])
    kw_plot_azav_grid.ncol = qgroup['ncol']
else:
    kw_plot_azav_grid.titles = parse_quantities(kw.qvals)[1]
    kw.groupname = input("choose a groupname to save your plot\n to not save it, enter 'nosave': ")

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# get data
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], 'AZ_Avgs')

print ('Getting data from ' + kw.the_file)
di = get_dict(kw.the_file)
vals = di['vals']
lut = di['lut']

# collect terms to plot
print ("plotting the following quantities:")
print ("qvals = ", kw.qvals)

terms = []
for qval in kw.qvals:
    the_term = get_term(dirname, vals, lut, qval, verbose=True)
    terms.append(the_term)

# for teq, make some modifications to the terms
if kw.groupname == 'teq':
    if advect_reference_state:
        terms.insert(3, terms[1]-terms[2]) # sum of cond + ref
        kw_plot_azav_grid.titles.insert(3, 'cond + ref')
        kw_plot_azav_grid.totsig.insert(3, 0)
        kw_plot_azav_grid.ncol += 1

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
ntheta = np.shape(vals)[0]
di_grid = get_grid_info(dirname, ntheta=ntheta)
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
