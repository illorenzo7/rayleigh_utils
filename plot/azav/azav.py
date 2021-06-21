##################################################################
# Routine to plot arbitrary variables (AZ_Avgs)
# Author: Loren Matilsky
# Created: 01/28/2019
##################################################################
# This script plots the quantitities specified by --qvals
# (default set in common)
##################################################################

import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
from azav_util import plot_azav_grid
from common import *
from plotcommon import *
from cla_util import *

# Read command-line arguments (CLAs)
args = sys.argv
if not '--qvals' in args:
    args += ['--qvals', 'v'] # make default qvals = v
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# get the data type we want; not all will be an AZ_Avgs file
dataname_list = dict({})
for ext in ['tot', 'pmp', 'ppm', 'mmm', 'mpp', 'ppp']:
    for ext2 in ['', 'r', 't', 'p']:
        dataname_list['meprodnum' + ext + ext2] = 'me_prod'
        dataname_list['meprodshear' + ext + ext2] = 'me_prod_shear'
for direc in ['r', 't', 'p']:
    dataname_list['ind' + direc + 'alt'] = 'induct_alt'
    dataname_list['ind' + direc + 'altnum'] = 'induct_alt_num'
for ext in ['mm', 'ms']:
    dataname_list['magtorque' + ext] = 'mag_torque'
for ext in ['tot', 'mmm', 'mpp']:
    dataname_list['meprodmean' + ext] = 'me_prod_mean'
dataname_list['ferraro'] = 'ferraro'

# See if magnetism is "on"
magnetism = clas0['magnetism']

# get desired quantities
qvals = clas['qvals']
titles = clas['titles']
ncol = clas['ncol']

if not clas0['tag'] == '': # the "tag" represents a quantity group
    clas0['tag'] = clas0['tag'][1:] # remove the prepending _

if clas0['tag'] in dataname_list.keys():
    dataname = dataname_list[clas0['tag']]
else:
    dataname = 'AZ_Avgs'

print ("plotting the following quantities:")
print ("qvals = " + arr_to_str(qvals, "%i"))
# get data
if 'the_file' in clas: 
    the_file = clas['the_file']
else:
    the_file = get_widest_range_file(clas0['datadir'], dataname)

print ('Getting quantities from ' + the_file)
di = get_dict(the_file)
vals = di['vals']
if dataname == 'AZ_Avgs':
    lut = di['lut']

# see if the user wants a separate plot of lat. averaged quantities
if 'shav' in clas:
    shav = True
else:
    shav = False

terms = []
for qval in qvals:
    if dataname == 'AZ_Avgs':
        terms.append(vals[:, :, lut[qval]])
    else:
        terms.append(vals[:, :, qval])

# make the main title
iter1, iter2 = get_iters_from_file(the_file)
time_string = get_time_string(dirname, iter1, iter2)
if not clas0['tag'] == '':
    mainlabel = clas0['tag']
else:
    mainlabel = 'quantities' 
maintitle = dirname_stripped + '\n' +\
        mainlabel + ' (zonally averaged)' + '\n' +\
        time_string

# Generate the figure using standard routine
di_grid = get_grid_info(dirname)
ncheby, domain_bounds = get_domain_bounds(dirname)
figs = plot_azav_grid (terms, di_grid['rr'], di_grid['cost'], maintitle=maintitle, tw=di_grid['tw'], domain_bounds=domain_bounds/rsun, **clas)
if shav:
    fig, av_fig = figs
else:
    fig = figs

# save the figure if tag was specified
if clas0['tag'] in ['v', 'b']: # these go in main directory
    plotdir = my_mkdir(clas0['plotdir'])
else:
    plotdir = my_mkdir(clas0['plotdir'] + 'azav/')
if not clas0['tag'] == '':
    savefile = plotdir + clas0['tag'] + '_azav-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
    print ('saving figure at ' + savefile)
    fig.savefig(savefile, dpi=300)

    if shav:
        av_savefile = plotdir + clas0['tag'] + '_shav-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
        print ('saving lat. avg. figure at ' + av_savefile)
        av_fig.savefig(av_savefile, dpi=300)

if clas0['showplot']:
    plt.show()

plt.close(fig)
if shav:
    plt.close(av_fig)
