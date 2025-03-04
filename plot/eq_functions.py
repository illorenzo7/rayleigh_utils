# Created: 05/03/2019
# Author: Loren Matilsky

import numpy as np
import matplotlib.pyplot as plt

import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
sys.path.append(os.environ['rapp'])

from common import *
from plotcommon import *
from cla_util import *
from reference_tools import equation_coefficients

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0.dirname
dirname_stripped = strip_dirname(dirname)

# allowed args + defaults
kw_default = dotdict()
kw_default.fname = None
kw_make_figure_default['margin_top_inches'] = 1.25
kw_default.update(kw_make_figure_default)
kw_default.update(kw_lineplot_default)

# change kw with clas
kw = update_dict(kw_default, clas)
kw_make_figure = update_dict(kw_make_figure_default, clas)
kw_lineplot = update_dict(kw_lineplot_default, clas)

# find bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# read equation_coefficients file
if kw.fname is None:
    kw.fname = 'equation_coefficients'
eq = equation_coefficients()
eq.read(dirname + '/' + kw.fname)

# things to plot and ylabels
f_dict = reverse_dict(eq.f_dict)
profiles = []
ylabels = []
for i in range(eq.nfunc):
    profiles.append(eq.functions[i, :])
    ylabels.append('f_%i = %s' %(i+1, f_dict[i+1]))

# create the plot; start with plotting all the energy fluxes
nplots = len(profiles)
ncol = 4
kw_make_figure.nplots = nplots
kw_make_figure.ncol = ncol
fig, axs, fpar = make_figure(**kw_make_figure)

# x label
kw_lineplot.xlabel = 'radius'

#for iplot in [0]:
for iplot in range(nplots):
    ax = axs.flatten()[iplot]
    kw_lineplot.ylabel = ylabels[iplot]
    lineplot(eq.radius, profiles[iplot], ax, **kw_lineplot)

# make title 
the_title = dirname_stripped + '\nequation_coefficients\n' +\
        ('nconst = %i\n' %eq.nconst) +\
        ('nfunc = %i\n' %eq.nfunc) +\
        ('version = %i' %eq.version)
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
