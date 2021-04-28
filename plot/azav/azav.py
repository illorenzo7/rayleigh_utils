##################################################################
# Routine to plot arbitrary variables (AZ_Avgs)
# Author: Loren Matilsky
# Created: 01/28/2019
##################################################################
# This script plots the quantitities specified by --qvals
# (default set in common)
##################################################################

import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
from azav_util import plot_azav_grid
from common import *
from plotcommon import *
from cla_util import *

# Read command-line arguments (CLAs)
clas = read_clas(sys.argv)
dirname = clas['dirname']
dirname_stripped = strip_dirname(dirname)

# get the data type we want; not all will be an AZ_Avgs file
dataname_list = dict({})
for ext in ['tot', 'pmp', 'ppm', 'mmm', 'mpp', 'ppp']:
    dataname_list['meprodnum' + ext] = 'me_prod'
    dataname_list['meprodshear' + ext] = 'me_prod_shear'
for direc in ['r', 't', 'p']:
    dataname_list['ind' + direc + 'alt'] = 'induct_alt'
    dataname_list['ind' + direc + 'altnum'] = 'induct_alt_num'
for ext in ['mm', 'ms']:
    dataname_list['magtorque' + ext] = 'mag_torque'
for ext in ['tot', 'mmm', 'mpp']:
    dataname_list['meprodmean' + ext] = 'me_prod_mean'

# See if magnetism is "on"
magnetism = get_parameter(dirname, 'magnetism')

# get desired quantities
qvals = clas['qvals']
if qvals is None:
    the_qgroup = get_quantity_group('v', magnetism)
    qvals = the_qgroup['qvals']
    titles = the_qgroup['titles']
    units = the_qgroup['units']
    ncol = the_qgroup['ncol']
    clas['tag'] = '_v' 
else:
    titles = clas['titles']
    units = clas['units']

if not clas['tag'] == '': # the "tag" represents a quantity group
    clas['tag'] = clas['tag'][1:] # remove the prepending _

if clas['tag'] in dataname_list.keys():
    dataname = dataname_list[clas['tag']]
else:
    dataname = 'AZ_Avgs'

print ("plotting the following quantities:")
print ("qvals = " + arr_to_str(qvals, "%i"))
the_file = clas['the_file']
if the_file is None:
    the_file = get_widest_range_file(clas['datadir'], dataname)
print ('Getting quantities from ' + the_file)
di = get_dict(the_file)
vals = di['vals']
if dataname == 'AZ_Avgs':
    lut = di['lut']

# see if the user wants a separate plot of lat. averaged quantities
latav = read_cla_arbitrary(sys.argv, 'latav', False)

terms = []
for qval in qvals:
    if dataname == 'AZ_Avgs':
        terms.append(vals[:, :, lut[qval]])
    else:
        terms.append(vals[:, :, qval])

# make the main title
iter1, iter2 = get_iters_from_file(the_file)
time_string = get_time_info(dirname, iter1, iter2)
if not clas['tag'] == '':
    mainlabel = clas['tag']
else:
    mainlabel = 'quantities' 
maintitle = dirname_stripped + '\n' +\
        mainlabel + ' (zonally averaged)' + '\n' +\
        time_string

# Generate the figure using standard routine
di_grid = get_grid_info(dirname)
figs = plot_azav_grid (terms, di_grid['rr'], di_grid['cost'], units=units, maintitle=maintitle, titles=titles,\
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
    plotboundary=clas['plotboundary'],\
    ncol=clas['ncol'],\
    fig_width_inches=clas['fig_width_inches'],\
    sub_width_inches=clas['sub_width_inches'],\
    latav=latav, tw=di_grid['tw'],\
    totsig=clas['totsig'])
if latav:
    fig, av_fig = figs
else:
    fig = figs

# save the figure if tag was specified
plotdir = my_mkdir(clas['plotdir'] + '/azav/')
if not clas['tag'] == '':
    savefile = plotdir + clas['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
    print ('saving figure at ' + savefile)
    fig.savefig(savefile, dpi=300)

    if latav:
        av_savefile = plotdir + clas['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '-latav' + '.png'
        print ('saving lat. avg. figure at ' + av_savefile)
        av_fig.savefig(av_savefile, dpi=300)

if clas['showplot']:
    plt.show()

plt.close(fig)
if latav:
    plt.close(av_fig)
