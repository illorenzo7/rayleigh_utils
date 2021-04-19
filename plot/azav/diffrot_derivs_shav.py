# Author: Loren Matilsky
# Created: 03/01/2021
# This script plots mean-shear induction terms in the meridional plane 
# breaks up into latitudinal and radial shear
# ...for the Rayleigh run directory indicated by [dirname]. 
# To use an AZ_Avgs file
# different than the one associated with the longest averaging range, use
# -usefile [complete name of desired AZ_Avgs file]
# Saves plot in
# [dirname]_induction_phi_mean_[first iter]_[last iter].png

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
from azav_util import plot_azav
from common import *

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'

# Read command-line arguments (CLAs)
nadd = 0
nsubset = []
torques_to_add = []
showplot = True
saveplot = True
plotcontours = True
plotlatlines = True
plotboundary = True
minmax = None
linthresh = None
linscale = None
minmaxrz = None
linthreshrz = None
linscalerz = None
the_file = get_widest_range_file(datadir, 'AZ_Avgs')
forced = False
rvals = []
rbcz = None
symlog = False
tag = ''

plotdir = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-minmaxrz':
        minmaxrz = float(args[i+1]), float(args[i+2])
    elif arg == '-rbcz':
        rbcz = float(args[i+1])
    elif arg == '-noshow':
        showplot = False
    elif arg == '-nosave':
        saveplot = False
    elif arg == '-usefile':
        the_file = args[i+1]
        the_file = the_file.split('/')[-1]
    elif arg == '-depths':
        strings = args[i+1].split()
        for st in strings:
            rval = ro - float(st)*d
            rvals.append(rval)
    elif arg == '-depthscz':
        rm = domain_bounds[1]
        dcz = ro - rm
        strings = args[i+1].split()
        for st in strings:
            rval = ro - float(st)*dcz
            rvals.append(rval)
    elif arg == '-depthsrz':
        rm = domain_bounds[1]
        drz = rm - ri
        strings = args[i+1].split()
        for st in strings:
            rval = rm - float(st)*drz
            rvals.append(rval)
    elif arg == '-rvals':
        rvals = []
        strings = args[i+1].split()
        for st in strings:
            rval = float(st)*rsun
            rvals.append(rval)
    elif arg == '-rvalscm':
        rvals = []
        strings = args[i+1].split()
        for st in strings:
            rval = float(st)
            rvals.append(rval)
    elif arg == '-symlog':
        symlog = True
    elif arg == '-linthresh':
        linthresh = float(args[i+1])
    elif arg == '-linscale':
        linscale = float(args[i+1])
    elif arg == '-linthreshrz':
        linthreshrz = float(args[i+1])
    elif arg == '-linscalerz':
        linscalerz = float(args[i+1])
    elif arg == '-nocontour':
        plotcontours = False
    elif arg == '-nobound':
        plotboundary = False
    elif arg == '-nolat':
        plotlatlines = False
    elif arg == '-add':
        loc_list = args[i+1].split()
        nadd += 1
        for j in range(len(loc_list)):
            loc_list[j] = int(loc_list[j])
        torques_to_add += loc_list
        nsubset.append(len(loc_list))
    elif arg == '-tag':
        tag = '_' + args[i+1]

# Get the terms:
print ('Getting terms from ' + datadir + the_file)
di = get_dict(datadir + the_file)

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']

# Get the time range in sec
t1 = translate_times(iter1, dirname, translate_from='iter')['val_sec']
t2 = translate_times(iter2, dirname, translate_from='iter')['val_sec']

# Get the baseline time unit
#rotation = get_parameter(dirname, 'rotation')
time_unit = compute_Prot(dirname)
Om0 = 2*np.pi/time_unit
time_label = r'$\rm{P_{rot}}$'

if plotdir is None:
    plotdir = dirname + '/plots/'
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)

# Get necessary grid info
rr = di['rr']
rr_2d = di['rr_2d']
cost = di['cost']
sint = di['sint']
cost_2d = di['cost_2d']
sint_2d = di['sint_2d']
tt_lat = di['tt_lat'] 
xx = di['xx'] 
nr, nt = di['nr'], di['nt'] 

# Get <v_phi>
vp = vals[:, :, lut[3]]

# get angular momentum viscous fluxes, and thus derivatives of <v_phi>
amom_visc_r = vals[:, :, lut[1813]]
amom_visc_t = vals[:, :, lut[1814]]

eq = get_eq(dirname)
nu = eq.nu.reshape((1, nr))
rho = eq.rho.reshape((1, nr))
mu = rho*nu
prefactor = -1./(rho*nu*rr_2d**2.*sint_2d**2.)

# get diffrot and its derivs, spherically averaged
gi = GridInfo(dirname + '/grid_info', '')
# make lat cutoff
lat_cutoff = 60.
ilat1 = np.argmin(np.abs(tt_lat + lat_cutoff))
ilat2 = np.argmin(np.abs(tt_lat - lat_cutoff))
nt_cut = ilat2 - ilat1 + 1
tw = gi.tweights[ilat1:ilat2+1]
tw /= np.sum(tw)
tw = tw.reshape((nt_cut, 1))
Om = vp/(rr_2d*sint_2d)
dOmdr_az = prefactor*amom_visc_r
dOmdt_az = prefactor*amom_visc_t
dOmdr = np.sqrt(np.sum(tw*(dOmdr_az**2.)[ilat1:ilat2+1,:], axis=0))
dOmdt = np.sqrt(np.sum(tw*(dOmdt_az**2.)[ilat1:ilat2+1,:], axis=0))
dOm =  np.sqrt(np.sum(tw*(dOmdr_az**2. + dOmdt_az**2.)[ilat1:ilat2+1,:], axis=0))

# get <v_phi> derivs
dvpdr_az = 1./rr_2d*(vp - amom_visc_r/mu/sint_2d)
dvpdt_az = 1./rr_2d/sint_2d*(cost_2d*vp - amom_visc_t/mu)
dvpdr = np.sqrt(np.sum(tw*dvpdr_az[ilat1:ilat2+1,:]**2., axis=0))
dvpdt = np.sqrt(np.sum(tw*dvpdt_az[ilat1:ilat2+1,:]**2., axis=0))
dvp =  np.sqrt(np.sum(tw*(dvpdr_az[ilat1:ilat2+1,:]**2. + dvpdt_az[ilat1:ilat2+1,:]**2.), axis=0))

terms = [dOmdr/(Om0/rsun), dOmdt/(Om0/rsun), dOm/(Om0/rsun), dOmdr/dOmdt,\
        dvpdr/Om0, dvpdt/Om0, dvp/Om0, dvpdr/dvpdt]
titles = ['dOm/dr', 'dOm/dT', 'magnitude/(Om0/rsun)', 'r/T ratio',\
    'dvp/dr', 'dvp/dT', 'magnitude/Om0', 'r/T ratio']

nplots = len(terms)
ncol = 4
nrow = np.int(np.ceil(nplots/ncol))

fig, axs = plt.subplots(nrow, ncol, figsize=(12., 3.*nrow))
iplot = 0
for irow in range(nrow):
    for icol in range(ncol):
        ax = axs[irow, icol]
        term = terms[iplot]
        ax.plot(rr, term)
        ax.set_title(titles[iplot])
        iplot += 1
savefile = plotdir + dirname_stripped + '_diffrot_derivs_shav_' +\
        str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'

plt.savefig(savefile, dpi=300)
plt.show()
