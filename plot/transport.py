# Author: Loren Matilsky
# Created: 12/19/2022
#
# Description: Script to plot radial diffusivity profiles 
# from the equation_coefficients file

import numpy as np
import matplotlib.pyplot as plt

import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])

from common import *
from plotcommon import *
from cla_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0.dirname
dirname_stripped = strip_dirname(dirname)
magnetism = clas0.magnetism

# allowed args + defaults
kw_default = dotdict()
kw_default.fname = None
kw_default.update(kw_make_figure_default)
kw_default.update(kw_lineplot_default)

# change kw with clas
kw = update_dict(kw_default, clas)
kw_make_figure = update_dict(kw_make_figure_default, clas)
kw_lineplot = update_dict(kw_lineplot_default, clas)

# find bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# read reference state
eq = get_eq(dirname, kw.fname, verbose=True)
print ("plotting transport-coefficient profiles")

# things to plot and ylabels
nplots = 4
nrow = 2
profiles = [eq.nu, eq.kappa, eq.dlnu, eq.dlnkappa]
ylabels = ['viscosity (' + r'$\nu$' + ')', 'thermometric conductivity (' + r'$\kappa$' + ')', r'$dln\nu/dr$', r'$dln\kappa/dr$']

# things change a bit with magnetism
if magnetism:
    nplots += 2
    profiles.insert(2, eq.eta)
    ylabels.insert(2, 'magnetic diffusivity (' + r'$\eta$' + ')')
    profiles.insert(5, eq.dlneta)
    ylabels.insert(5, r'$dln\eta/dr$')

# Create the plot; start with plotting all the energy fluxes
kw_make_figure.nplots = nplots
kw_make_figure.nrow = nrow
fig, axs, fpar = make_figure(**kw_make_figure)

# x label
kw_lineplot.xlabel = 'radius (r)'

#for iplot in [0]:
for iplot in range(nplots):
    ax = axs.flatten()[iplot]
    kw_lineplot.ylabel = ylabels[iplot]
    lineplot(eq.rr, profiles[iplot], ax, **kw_lineplot)

# make title 
the_title = dirname_stripped + '\ntransport coefficients (diffusivities)'
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
