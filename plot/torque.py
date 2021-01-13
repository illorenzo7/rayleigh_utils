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
from azav_util import plot_azav
from common import *
from read_inner_vp import read_inner_vp
from read_eq_vp import read_eq_vp

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

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
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
forced = False
rvals = []
rbcz = None
symlog = False
tag = ''

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
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
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
    elif arg == '-forced':
        forced = True
    elif arg == '-depths':
        rvals = []
        strings = args[i+1].split()
        for st in strings:
            rval = ro - float(st)*d
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

# See if magnetism is "on"
try:
    magnetism = get_parameter(dirname, 'magnetism')
except:
    magnetism = False # if magnetism wasn't specified, it must be "off"

# Get the torques:
print ('Getting torques from ' + datadir + AZ_Avgs_file)
di = get_dict(datadir + AZ_Avgs_file)

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']

# Get the time range in sec
t1 = translate_times(iter1, dirname, translate_from='iter')['val_sec']
t2 = translate_times(iter2, dirname, translate_from='iter')['val_sec']

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

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

max_sig = max(np.std(torque_rs), np.std(torque_mc), np.std(torque_visc))

if magnetism:
    ind_Maxwell_mean = lut[1805]
    ind_Maxwell_rs = lut[1806]
    
    torque_Maxwell_mean = vals[:, :, ind_Maxwell_mean]
    torque_Maxwell_rs = vals[:, :, ind_Maxwell_rs]
    
    torque_tot += torque_Maxwell_mean + torque_Maxwell_rs
    
    max_sig = max(max_sig, np.std(torque_Maxwell_mean),\
                  np.std(torque_Maxwell_rs))

if forced: # compute torque associated with forcing function
    mean_vp = vals[:, :, lut[3]]
    tacho_r = get_parameter(dirname, 'tacho_r')
    print ("read tacho_r = %1.2e" %tacho_r)
    tacho_dr = get_parameter(dirname, 'tacho_dr')
    tacho_tau = get_parameter(dirname, 'tacho_tau')
    forcing = np.zeros((nt, nr))
    eq = get_eq(dirname)
    rho = eq.rho

    if os.path.exists(dirname + '/inner_vp'):
        print ("inner_vp file exists, so I assume you have a forcing function which\n quartically matches on to a CZ differential rotation")
        inner_vp = read_inner_vp(dirname + '/inner_vp')
        for it in range(nt):
            for ir in range(nr):
                if rr[ir] <= tacho_r - tacho_dr*rr[0]:
                    desired_vp = 0.0
                elif rr[ir] <= tacho_r:
                    desired_vp = inner_vp[it]*(1.0 - ( (rr[ir] - tacho_r)/(tacho_dr*rr[0]) )**2)**2
                else:
                    desired_vp = mean_vp[it, ir] # No forcing in CZ
                forcing[it, ir] = -rho[ir]*(mean_vp[it, ir] - desired_vp)/tacho_tau
    elif os.path.exists(dirname + '/eq_vp'):
        print ("eq_vp file exists, so I assume you have a forcing function which\n quartically matches on to a CZ differential rotation\n with viscous-torque-free buffer zone")
        tacho_r2 = get_parameter(dirname, 'tacho_r2')
        i_tacho_r = np.argmin(np.abs(rr - tacho_r))
        print ("read tacho_r2 = %1.2e" %tacho_r2)
        eq_vp = read_eq_vp(dirname + '/eq_vp', nt, nr)
        for it in range(nt):
            for ir in range(nr):
                if rr[ir] <= tacho_r - tacho_dr*rr[0]:
                    desired_vp = 0.0
                elif rr[ir] <= tacho_r:
                    desired_vp = eq_vp[it, i_tacho_r]*(1.0 - ( (rr[ir] - tacho_r)/(tacho_dr*rr[0]) )**2)**2
                elif rr[ir] <= tacho_r2:
                    desired_vp = eq_vp[it, ir]
                else:
                    desired_vp = mean_vp[it, ir] # No forcing in CZ
                forcing[it, ir] = -rho[ir]*(mean_vp[it, ir] - desired_vp)/tacho_tau
    else:
        forcing_coeff = -rho/tacho_tau*0.5*(1.0 - np.tanh((rr - tacho_r)/(tacho_dr*rr[0])))
        forcing = forcing_coeff.reshape((1, nr))*mean_vp

    # convert forcing function into a torque
    torque_forcing = di['xx']*forcing
    torque_tot += torque_forcing
    
    max_sig = max(max_sig, np.std(torque_forcing))

# Set up the actual figure from scratch
fig_width_inches = 7. # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1./8. # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)
margin_top_inches = 1 # wider top margin to accommodate subplot titles AND metadata
margin_subplot_top_inches = 1/4 # margin to accommodate just subplot titles
nplots = 4 + 2*magnetism + 1*forced + nadd
ncol = 3 # put three plots per row
nrow = np.int(np.ceil(nplots/3))

subplot_width_inches = (fig_width_inches - (ncol + 1)*margin_inches)/ncol
    # Make the subplot width so that ncol subplots fit together side-by-side
    # with margins in between them and at the left and right.
subplot_height_inches = 2*subplot_width_inches # Each subplot should have an
    # aspect ratio of y/x = 2/1 to accommodate meridional planes. 
fig_height_inches = margin_top_inches + nrow*(subplot_height_inches +\
        margin_subplot_top_inches + margin_bottom_inches)
fig_aspect = fig_height_inches/fig_width_inches

# "Margin" in "figure units"; figure units extend from 0 to 1 in BOTH 
# directions, so unitless dimensions of margin will be different in x and y
# to force an equal physical margin
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches
margin_subplot_top = margin_subplot_top_inches/fig_height_inches

# Subplot dimensions in figure units
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

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
    ind_insert += 2
if forced:
    torques.insert(ind_insert, torque_forcing)
    titles.insert(ind_insert, r'$\tau_{\rm{forcing}}$')
    ind_insert =+ 1
    
if nadd > 0:
    jstart = 0
    for i in range(nadd):
        ind = torques_to_add[jstart]
        torque_sum = np.copy(torques[ind])
        title_sum = titles[ind]
        for j in range(jstart + 1, jstart + nsubset[i]):
            ind = torques_to_add[j]
            torque_sum += torques[ind]
            title_sum += r'$\ +\ $' + titles[ind]
        jstart += nsubset[i]

        torques.insert(ind_insert, torque_sum)
        titles.insert(ind_insert, title_sum)
        ind_insert += 1

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

for iplot in range(nplots):
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1 - margin_top - subplot_height - margin_subplot_top -\
            (iplot//ncol)*(subplot_height + margin_subplot_top +\
            margin_bottom)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    plot_azav (torques[iplot], rr, cost, fig=fig, ax=ax, units=units,\
           minmax=minmax, plotcontours=plotcontours, rvals=rvals,\
           minmaxrz=minmaxrz, rbcz=rbcz, symlog=symlog,\
    linthresh=linthresh, linscale=linscale, linthreshrz=linthreshrz,\
    linscalerz=linscalerz, plotlatlines=plotlatlines, plotboundary=plotboundary)

    ax.set_title(titles[iplot], verticalalignment='bottom', **csfont)

# Label averaging interval
if rotation:
    time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.1f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'
else:
    time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.3f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'

# Put some metadata in upper left
fsize = 12
fig.text(margin_x, 1 - 0.1*margin_top, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.3*margin_top, 'Torque balance (zonally averaged)',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.5*margin_top, time_string,\
         ha='left', va='top', fontsize=fsize, **csfont)

savefile = plotdir + dirname_stripped + '_torque_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + tag + '.png'

if saveplot:
    print ('Saving torques at ' + savefile)
    plt.savefig(savefile, dpi=300)
if showplot:
    plt.show()
plt.close()
