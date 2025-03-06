# Author: Loren Matilsky
# Created: 12/19/2022
#
# Description: Script to plot varius length scales as functions of
# radius using Shell_Avgs data

import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
from common import *
from plotcommon import *
from cla_util import *

# get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# allowed args + defaults
kw_default = dict({'the_file': None})
kw_default.update(kw_make_figure_default)
#kw_lineplot_default['legfrac'] = 0.3
kw_lineplot_default['plotleg'] = True
kw_lineplot_default['ncolleg'] = 2
kw_lineplot_default['log'] = True
kw_default.update(kw_lineplot_default)
kw = update_dict(kw_default, clas)
kw_make_figure = update_dict(kw_make_figure_default, clas)
kw_lineplot = update_dict(kw_lineplot_default, clas)

di_ls = length_scales(dirname)
rr = di_ls.rr
kw_lineplot.labels = ['|v|/|om|', '|<v>|/|<om>|', "|v'|/|om'|"]
profiles = [di_ls.v, di_ls.vmean, di_ls.vfluc]

magnetism = clas0.magnetism
if magnetism:
    profiles += [di_ls.b, di_ls.bmean, di_ls.bfluc]
    kw_lineplot.labels += ['|B|/|J|', '|<B>|/|<J>|', "|B'|/|J'|"]

# create the plot
fig, axs, fpar = make_figure(**kw_make_figure)
ax = axs[0,0]
kw_lineplot.xlabel = 'radius'
kw_lineplot.ylabel = 'length scale'

lineplot(rr, profiles, ax, **kw_lineplot)

# make title 
iter1, iter2 = di_ls.iter1, di_ls.iter2
time_string = get_time_string(dirname, iter1, iter2) 
the_title = dirname_stripped + '\n' +  'fluid length scales' + '\n' + time_string
if kw.mark_bcz:
    the_title += ('\n' + r'$r_{BCZ} = %1.3e$' %rbcz_est)
    the_title += ('\n' + r'$r_{os} = %1.3e$' %rov_est)

margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
fig.text(margin_x, 1 - margin_y, the_title,\
         ha='left', va='top', fontsize=default_titlesize)


if clas0['saveplot']:
    plotdir = my_mkdir(clas0['plotdir'])
    savefile = plotdir + clas0['routinename'] + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
