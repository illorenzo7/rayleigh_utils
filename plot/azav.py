##################################################################
# Routine to plot arbitrary variables (AZ_Avgs)
# Author: Loren Matilsky
# Created: 01/28/2019
##################################################################
# This script plots the quantitities specified by --qvals
# (default set in common)
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
magnetism = get_parameter(dirname, 'magnetism')

# get desired quantities
qvals = clas['qvals']
titles = clas['titles']
units = clas['units']
print ("plotting the following quantities:")
print ("qvals = " + arr_to_str(qvals, "%i"))
the_file = clas['the_file']
if the_file is None:
    the_file = get_widest_range_file(datadir, 'AZ_Avgs')
print ('Getting quantities from ' + datadir + the_file)
di = get_dict(datadir + the_file)

vals = di['vals']
lut = di['lut']

plotdir = clas['plotdir']
if plotdir is None:
    plotdir = dirname + '/plots/'
make_plotdir(plotdir)

terms = []
for qval in qvals:
    terms.append(vals[:, :, lut[qval]])

# make the main title
iter1, iter2 = get_iters_from_file(the_file)
time_string = get_time_info(dirname, iter1, iter2)
maintitle = dirname_stripped + '\n' +\
        'Quantities (zonally averaged)' + '\n' +\
        time_string

# Generate the figure using standard routine
di_grid = get_grid_info(dirname)
fig = plot_azav_grid (terms, di_grid['rr'], di_grid['cost'], units=units, maintitle=maintitle, titles=titles,\
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
    ncol=clas['ncol'])

# save the figure if tag was specified
tag = clas['tag']
if not tag == '':
    tag = tag[1:] # (in this case, remove the prepending underscore
    savefile = plotdir + tag + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas['showplot']:
    plt.show()
plt.close()
