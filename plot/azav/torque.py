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
import pickle
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
from azav_util import plot_azav_grid
from common import *
from cla_util import *

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
#already exist    
datadir = dirname + '/data/'

# Read command-line arguments (CLAs)
args = sys.argv[2:]
clas = read_clas(dirname, args)

# See if magnetism is "on"
try:
    magnetism = get_parameter(dirname, 'magnetism')
except:
    magnetism = False # if magnetism wasn't specified, it must be "off"

# Get the torques
the_file = clas['the_file']
if the_file is None:
    the_file = get_widest_range_file(datadir, 'AZ_Avgs')
print ('Getting torques from ' + datadir + the_file)
di = get_dict(datadir + the_file)

vals = di['vals']
lut = di['lut']

plotdir = make_plotdir(dirname, clas['plotdir'], '/plots/azav/')

torque_rs, torque_mc, torque_visc = -vals[:, :, lut[1801]],\
        -vals[:, :, lut[1802]] + vals[:, :, lut[1803]],\
        vals[:, :, lut[1804]]
torque_tot = torque_rs + torque_mc + torque_visc

if magnetism:
    torque_Maxwell_mean = vals[:, :, lut[1805]]
    torque_Maxwell_rs = vals[:, :, lut[1806]]
    torque_tot += torque_Maxwell_mean + torque_Maxwell_rs

torques = [torque_rs, torque_mc, torque_visc, torque_tot]
titles = [r'$\tau_{\rm{rs}}$', r'$\tau_{\rm{mc}}$', r'$\tau_{\rm{v}}$',\
          r'$\tau_{\rm{tot}}$']
units = r'$\rm{g}\ \rm{cm}^{-1}\ \rm{s}^{-2}$'

ind_insert = 3
if magnetism:
    torques.insert(ind_insert, torque_Maxwell_mean)
    torques.insert(ind_insert + 1, torque_Maxwell_rs)
    titles.insert(ind_insert, r'$\tau_{\rm{mm}}$') 
    titles.insert(ind_insert + 1, r'$\tau_{\rm{ms}}$')
    #ind_insert += 2

# make the main title
iter1, iter2 = get_iters_from_file(the_file)
time_string = get_time_info(dirname, iter1, iter2)
maintitle = dirname_stripped + '\n' +\
        'Torque balance (zonally averaged)' + '\n' +\
        time_string

# Generate the figure using standard routine
di_grid = get_grid_info(dirname)
fig = plot_azav_grid (torques, di_grid['rr'], di_grid['cost'], units=units, maintitle=maintitle, titles=titles,\
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
savefile = plotdir + 'torque' + clas['tag'] + '-' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + '.png'

if clas['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas['showplot']:
    plt.show()
plt.close()
