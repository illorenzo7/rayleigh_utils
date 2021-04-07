# Author: Loren Matilsky
# Created: 01/28/2019
# This script plots the axial torques in the meridional plane (viscous, 
# Meridional Circ., Reynolds stress, and Maxwell torques (mean and turbulent) 
# if applicablein the meridional plane 
# ...for the Rayleigh run directory indicated by [dirname]. To use an AZ_Avgs file
# different than the one associated with the longest averaging range, use
# -usefile [complete name of desired AZ_Avgs file]
# Saves plot in
# [dirname]_torque_[first iter]_[last iter].png

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
from read_inner_vp import read_inner_vp
from read_eq_vp import read_eq_vp

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'

# Read command-line arguments (CLAs)
args = sys.argv[2:]
clas = read_clas(args)
rvals = read_rvals(dirname, args)
#nargs = len(args)
#for i in range(nargs):
#    arg = args[i]
#    if arg == '-plotdir':
#        plotdir = args[i+1]

# See if magnetism is "on"
try:
    magnetism = get_parameter(dirname, 'magnetism')
except:
    magnetism = False # if magnetism wasn't specified, it must be "off"

# Get the torques:
thefile = clas['thefile'].val
if thefile is None:
    thefile = get_widest_range_file(datadir, 'AZ_Avgs')
print ('Getting torques from ' + datadir + thefile)
di = get_dict(datadir + thefile)

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']

plotdir = clas['plotdir'].val
if plotdir is None:
    plotdir = dirname + '/plots/'
make_plotdir(plotdir)

# Get necessary grid info
rr = di['rr']
cost = di['cost']
sint = di['sint']
tt_lat = di['tt_lat']
xx = di['xx']
nr, nt = di['nr'], di['nt']

ind_pp = lut[1801]
ind_mm = lut[1802]
ind_cor = lut[1803]
ind_visc = lut[1804]

torque_rs, torque_mc, torque_visc = -vals[:, :, ind_pp],\
        -vals[:, :, ind_mm] + vals[:, :, ind_cor],\
        vals[:, :, ind_visc]
torque_tot = torque_rs + torque_mc + torque_visc

if magnetism:
    ind_Maxwell_mean = lut[1805]
    ind_Maxwell_rs = lut[1806]
    
    torque_Maxwell_mean = vals[:, :, ind_Maxwell_mean]
    torque_Maxwell_rs = vals[:, :, ind_Maxwell_rs]
    
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
time_string, time_unit, time_label, time_key =\
        get_time_info(dirname, iter1, iter2)
maintitle = dirname_stripped + '\n' +\
        'Torque balance (zonally averaged)' + '\n' +\
        time_string

# Generate the figure using standard routine
fig = plot_azav_grid (torques, rr, cost, units=units, maintitle=maintitle, titles=titles,\
        minmax=clas['minmax'].val,\
        plotcontours=clas['plotcontours'].val,\
        rvals=rvals,\
        minmaxrz=clas['minmaxrz'].val,\
        rbcz=clas['rbcz'].val,\
        symlog=clas['symlog'].val,\
    linthresh=clas['linthresh'].val,\
    linscale=clas['linscale'].val,\
    linthreshrz=clas['linthreshrz'].val,\
    linscalerz=clas['linscalerz'].val,\
    plotlatlines=clas['plotlatlines'].val,\
    plotboundary=clas['plotboundary'].val)

# save the figure
savefile = plotdir + dirname_stripped + '_torque_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + clas['tag'].val + '.png'

if clas['saveplot'].val:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas['showplot'].val:
    plt.show()
plt.close()
