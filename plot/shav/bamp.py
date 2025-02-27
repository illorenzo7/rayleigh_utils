# Author: Loren Matilsky
# Created: 10/31/2019

# Import relevant modules
import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.environ['rapl'])
sys.path.append(os.environ['raco'])
from common import *
from plotcommon import *
from cla_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# allowed args + defaults
kw_default = dict({'the_file': None, 'type': 'fluc'})
kw_default.update(kw_make_figure_default)
kw_lineplot_default['legfrac'] = 0.25
kw_lineplot_default['buff_ignore'] = buff_frac # ignore nastiness 
# at endpoints
kw_default.update(kw_lineplot_default)
kw = update_dict(kw_default, clas)
kw_make_figure = update_dict(kw_make_figure_default, clas)
kw_lineplot = update_dict(kw_lineplot_default, clas)
if not kw.xcut is None: # make room for label on right
    kw_make_figure.sub_margin_right_inches = default_margin_xlabel
    kw_make_figure.margin_top_inches += default_line_height
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# get data and grid
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], 'Shell_Avgs')
print ('Getting data from ' + kw.the_file)
di = get_dict(kw.the_file)
vals = di['vals']
lut = di['lut']
di_grid = get_grid_info(dirname)
rr = di_grid['rr']
nr = di_grid['nr']

# get c_4/2
eq = get_eq(dirname)
c4 = eq.constants[3]
coeff = c4/2.

# Convective velocity amplitudes, get these from me
rme_fluc = vals[:, 0, lut[1110]]
tme_fluc = vals[:, 0, lut[1111]]
pme_fluc = vals[:, 0, lut[1112]]

# Total velocity amplitudes, get these from me
rme = vals[:, 0, lut[1102]]
tme = vals[:, 0, lut[1103]]
pme = vals[:, 0, lut[1104]]

# mean velocity amplitudes from tot - fluc
rme_mean = rme - rme_fluc
tme_mean = tme - tme_fluc
pme_mean = pme - pme_fluc

if kw.type == 'fluc':
    rme_use = rme_fluc
    tme_use = tme_fluc
    pme_use = pme_fluc
elif kw.type == 'tot':
    rme_use = rme
    tme_use = tme
    pme_use = pme
elif kw.type == 'mean':
    rme_use = rme_mean
    tme_use = tme_mean
    pme_use = pme_mean
elif kw.type == 'rat': # compute a bunch of ratios
    bsq = (rme + tme + pme)/coeff # total
    bsq_fluc = (rme_fluc + tme_fluc + pme_fluc)/coeff # convect.
    bsq_r = rme/coeff # vertical total
    bsq_h = (tme + pme)/coeff # horizontal convect.
    bsq_tor = pme_mean/coeff # diff. rot.
    bsq_pol = (rme_mean + tme_mean)/coeff # mer. circ.

    amp_b = np.sqrt(bsq)
    amp_fluc = np.sqrt(bsq_fluc)
    amp_r = np.sqrt(bsq_r)
    amp_h = np.sqrt(bsq_h)
    amp_tor = np.sqrt(bsq_tor)
    amp_pol = np.sqrt(bsq_pol)
else:
    print("type must be fluc, mean, tot, or rat")
    print("exiting now")
    sys.exit()

# compute desired profiles
if kw.type == 'rat':
    profiles = [amp_fluc/amp_b, amp_tor/amp_b, amp_pol/amp_b, amp_r/amp_h]
    kw_lineplot.labels = [r'$|\mathbf{B}^\prime|/|\mathbf{B}|$', r'$|\langle\mathbf{B}_{tor}\rangle|/|\mathbf{B}|$', r'$|\langle\mathbf{B}_{pol}\rangle|/|\mathbf{B}|$', r'$|B_r|/|\mathbf{B}_h|$']
else:
    bsq_r = rme_use/coeff
    bsq_t = tme_use/coeff
    bsq_p = pme_use/coeff
    bsq = bsq_r + bsq_t + bsq_p

    amp_b = np.sqrt(bsq)
    amp_br = np.sqrt(bsq_r)
    amp_bt = np.sqrt(bsq_t)
    amp_bp = np.sqrt(bsq_p)

    profiles = [amp_br, amp_bt, amp_bp, amp_b]
    kw_lineplot.labels = [r'$B_r$', r'$B_\theta$', r'$B_\phi$', r'$\mathbf{B}$']

# Create the plot; start with plotting all the energy fluxes
fig, axs, fpar = make_figure(**kw_make_figure)
ax = axs[0,0]

# x and y labels
kw_lineplot.xlabel = 'radius'
kw_lineplot.ylabel = 'rms velocity'
kw_lineplot.plotleg = True
lineplot(rr, profiles, ax, **kw_lineplot)

# make title 
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2) 
the_title = dirname_stripped + '\n' +  ' velocity amplitudes (' +\
        kw.type + ')' + '\n' + time_string
ax.set_title(the_title, fontsize=default_titlesize)

# save the figure
plotdir = my_mkdir(clas0['plotdir'])
savefile = plotdir + clas0['routinename'] + '_' + kw.type + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
