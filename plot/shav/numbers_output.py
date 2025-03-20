# Created: 05/03/2019 # Author: Loren Matilsky

import numpy as np
import matplotlib.pyplot as plt

import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])

from common import *
from plotcommon import *
from fluid_numbers import *
from cla_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0.dirname
dirname_stripped = strip_dirname(dirname)
magnetism = clas0.magnetism
rotation = clas0.rotation

# allowed args + defaults
kw_default = dict({'the_file': None, 'the_file_az': None, 'sd': None})
kw_default.update(kw_make_figure_default)
kw_default.update(kw_lineplot_default)

# change kw with clas
kw = update_dict(kw_default, clas)
kw_make_figure = update_dict(kw_make_figure_default, clas)
kw_lineplot = update_dict(kw_lineplot_default, clas)

# find bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# deal with shell depth (by default use whole shell)
shell_depth = clas.sd
if shell_depth is None:
    rmin, rmax = interpret_rvals(dirname, ['rmin', 'rmax'])
    shell_depth = rmax - rmin

# get the output numbers
di = get_numbers_output(dirname, shell_depth, kw.the_file, kw.the_file_az)
rr = get_grid_info(dirname)['rr']

# number of plots (number of groups of output numbers)
nplots = numbers_output_ngroup
if rotation:
    nplots += numbers_output_ngroup_rot
if magnetism:
    nplots += numbers_output_ngroup_mag
ncol = 4
kw_make_figure.nplots = nplots
kw_make_figure.ncol = ncol
fig, axs, fpar = make_figure(**kw_make_figure)

# x label
kw_lineplot.xlabel = 'radius'

for iplot in range(nplots):
    # collect profiles between line breaks
    if iplot == 0:
        iprofmin = 0
    else:
        iprofmin = linebreaks_output[iplot-1]
    if iplot == nplots - 1:
        iprofmax = len(di)
    else:
        iprofmax = linebreaks_output[iplot]

    keys_loc = list(di.keys())[iprofmin:iprofmax]
    labels = []
    profiles = []

    for key in keys_loc:
        labels.append(numbers_output_def[key])
        profiles.append(di[key])

    # plot the profiles on same panel
    ax = axs.flatten()[iplot]
    kw_lineplot.ylabel = numbers_output_groups[iplot]
    kw_lineplot.labels = labels
    kw_lineplot.plotleg = True
    kw_lineplot.ncolleg = 1
    lineplot(rr, profiles, ax, **kw_lineplot)

# make title 
the_title = dirname_stripped + '\nNON-DIMENSIONAL NUMBERS (radial profiles)'
margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
fig.text(margin_x, 1 - margin_y, the_title,\
         ha='left', va='top', fontsize=default_titlesize)

# save the figure, maybe
if clas0['saveplot']:
    plotdir = my_mkdir(clas0['plotdir'])
    savefile = plotdir + clas0['routinename'] + clas0['tag'] + '.png'
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
