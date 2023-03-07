# Author: Loren Matilsky
# Created: 12/19/2022
#
# Description: Script to plot radial thermodynamic profiles 
# (the reference state) from the equation_coefficients file

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

# allowed args + defaults
kwargs_default = dotdict()
kwargs_default.fname = None
kwargs_default.update(make_figure_kwargs_default)
kwargs_default.update(lineplot_kwargs_default)

# change kwargs with clas
kw = update_dict(kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)
kw_lineplot = update_dict(lineplot_kwargs_default, clas)

# find bad keys
find_bad_keys(kwargs_default, clas, clas0['routinename'], justwarn=True)

# read reference state
eq = get_eq(dirname, kw.fname)

# things to plot and ylabels
profiles = [eq.grav, eq.dsdr, eq.heat,\
        eq.rho, eq.tmp, eq.dlnrho,\
        eq.d2lnrho, eq.dlntmp, eq.prs]
ylabels = ['gravity (g)', r'$d\overline{S}/dr$',  'heating (Q)',\
        'density (' + r'$\overline{\rho}$' + ')', 'temperature (' + r'$\overline{T}$' + ')', r'$dln\overline{\rho}/dr$', \
        r'$d^2ln\rho/dr^2$', r'$dln\overline{T}/dr$', 'pressure (' + r'$\overline{P}=\overline{\rho}\mathcal{R}\overline{T}$' + ')']
count = 0
if eq.reference_type in [2, 4]:
    profiles.insert(2, eq.nsq)
    ylabels.insert(2, r'$N^2=(g/c_p)d\overline{S}/dr$')
    count += 1
if not close_to_zero(eq.dlnrho):
    profiles.insert(6 + count, -1.0/eq.dlnrho)
    ylabels.insert(6 + count, r'$H_\rho=-(dln\overline{\rho}/dr)^{-1}$')

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
    lineplot(eq.rr, profiles[iplot], ax, **kw_lineplot)

# make title 
the_title = dirname_stripped + '\nbackground reference state' +\
        '\nreference_type = %i' %eq.reference_type
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
