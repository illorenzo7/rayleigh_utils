##################################################################
# Routine to plot torques (AZ_Avgs)
# Author: Loren Matilsky
# Created: 01/28/2019
##################################################################
# This script plots the axial torques in the meridional plane (viscous, 
# Meridional Circ., Reynolds stress, and possibly
# Maxwell torques (mean and turbulent) 
##################################################################

import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
#sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
from azav_util import plot_azav_grid
from varprops import utype
from common import *
from cla_util import *

# Get directory name and stripped_dirname for plotting purposes
args = sys.argv
clas = read_clas(args)
dirname = clas['dirname']
dirname_stripped = strip_dirname(dirname)

# See if magnetism is "on"
magnetism = get_parameter(dirname, 'magnetism')

# Get the torques
the_file = clas['the_file']
if the_file is None:
    the_file = get_widest_range_file(clas['datadir'], 'AZ_Avgs')
print ('Getting torques from ' + the_file)
di = get_dict(the_file)
vals = di['vals']
lut = di['lut']

torques = [-vals[:, :, lut[1801]], -vals[:, :, lut[1802]] + vals[:, :, lut[1803]], vals[:, :, lut[1804]]] # rs, mc, visc
titles = [r'$\tau_{\rm{rs}}$', r'$\tau_{\rm{mc}}$', r'$\tau_{\rm{v}}$']
units = utype['torque']

if magnetism:
    torques += [vals[:, :, lut[1805]], vals[:, :, lut[1806]]] # mm, ms
    titles += [r'$\tau_{\rm{mm}}$', r'$\tau_{\rm{ms}}$']

# get total torque
torque_tot = np.zeros_like(torques[0])
for torque in torques:
    torque_tot += torque
torques.append(torque_tot)
titles.append(r'$\tau_{\rm{tot}}$')

# make the main title
iter1, iter2 = get_iters_from_file(the_file)
time_string = get_time_info(dirname, iter1, iter2)
maintitle = dirname_stripped + '\n' +\
        'Torque balance (zonally averaged)' + '\n' +\
        time_string

# Generate the figure using standard routine
di_grid = get_grid_info(dirname)
fig = plot_azav_grid (torques, di_grid['rr'], di_grid['cost'],\
        units=units, maintitle=maintitle, titles=titles,\
        minmax=clas['minmax'],\
        plotcontours=clas['plotcontours'],\
        rvals=clas['rvals'],\
        minmaxrz=clas['minmaxrz'],\
        rbcz=clas['rbcz'],\
        symlog=clas['symlog'],\
    linthresh=clas['linthresh'],\
    linscale=clas['linscale'],\
    linthreshrz=clas['linthreshrz'],\
    linscalerz=clas['linscalerz'],\
    plotlatlines=clas['plotlatlines'],\
    plotboundary=clas['plotboundary'])

# save the figure
plotdir = my_mkdir(clas['plotdir'] + '/azav/')
savefile = plotdir + clas['routinename'] + clas['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas['showplot']:
    plt.show()
plt.close()
