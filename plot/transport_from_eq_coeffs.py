# Created: 05/03/2019
# Author: Loren Matilsky

import matplotlib as mpl
import numpy as np
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
#import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from get_parameter import get_parameter

from common import strip_dirname, rsun
from reference_tools import equation_coefficients

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# defaults
ylog = False

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-log':
        ylog = True

eq = equation_coefficients()
eq.read(dirname + '/equation_coefficients')
r = eq.radius
nu = eq.functions[2, :]
dlnu = eq.functions[10, :]
kappa = eq.functions[4, :]
dlnkappa = eq.functions[11, :]

magnetism = get_parameter(dirname, 'magnetism')
if magnetism:
    eta = eq.functions[6, :]
    dlneta = eq.functions[12, :]

if magnetism:
    figsize = 7., 10.
    fig, axs = plt.subplots(3, 2, figsize=figsize, sharex=True)
else:
    figsize = 7., 7.
    fig, axs = plt.subplots(2, 2, figsize=figsize, sharex=True)

# Plot transport coeff and log derivatives one-by-one
lw = 1.
axs[0,0].plot(r/rsun, nu, linewidth=lw)
#    axs[1,0].yaxis.set_major_formatter(yfmt)
axs[0,0].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
axs[0,0].set_ylabel(r'$\nu(r)\ $' +  r'$\rm{[cm^2\ s^{-1}]}$')
if ylog:
    axs[0,0].set_yscale('log')

axs[0,1].plot(r/rsun, dlnu, linewidth=lw)
#    axs[1,0].yaxis.set_major_formatter(yfmt)
axs[0,1].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
axs[0,1].set_ylabel(r'$d\ln\nu/dr\ $' +  r'$\rm{[cm^{-1}]}$')

axs[1,0].plot(r/rsun, kappa, linewidth=lw)
#    axs[1,0].yaxis.set_major_formatter(yfmt)
axs[1,0].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
axs[1,0].set_ylabel(r'$\kappa(r)\ $' +  r'$\rm{[cm^2\ s^{-1}]}$')
if ylog:
    axs[1,0].set_yscale('log')

axs[1,1].plot(r/rsun, dlnu, linewidth=lw)
#    axs[1,0].yaxis.set_major_formatter(yfmt)
axs[1,1].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
axs[1,1].set_ylabel(r'$d\ln\kappa/dr\ $' +  r'$\rm{[cm^{-1}]}$')

# Get ticks everywhere
plt.sca(axs[0,0])
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

plt.sca(axs[1,0])
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

plt.sca(axs[1,0])
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

plt.sca(axs[1,1])
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')    


if magnetism:
    axs[2,0].plot(r/rsun, eta, linewidth=lw)
    #    axs[1,0].yaxis.set_major_formatter(yfmt)
    axs[2,0].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[2,0].set_ylabel(r'$\nu(r)\ $' +  r'$\rm{[cm^2\ s^{-1}]}$')
    if ylog:
        axs[2,0].set_yscale('log')

    axs[2,1].plot(r/rsun, dlneta, linewidth=lw)
    #    axs[1,0].yaxis.set_major_formatter(yfmt)
    axs[2,1].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[2,1].set_ylabel(r'$d\ln\eta/dr\ $' +  r'$\rm{[cm^{-1}]}$')
    
    # ticks everywhere
    plt.sca(axs[2,0])
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    plt.sca(axs[2,1])
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')    

# Put some metadata in the title
axs[0,0].set_title(dirname_stripped +'\n' + 'transport coefficients',\
        ha='center', pad=20)

# Label the x axes of the bottom plots
xlabel = r'$r/R_\odot$'
if magnetism:
    axs[2,0].set_xlabel(xlabel)
    axs[2,1].set_xlabel(xlabel)
else:
    axs[1,0].set_xlabel(xlabel)
    axs[1,1].set_xlabel(xlabel)

# Set the x axis limits
axs[0,0].set_xlim(np.min(r)/rsun, np.max(r)/rsun)
plt.subplots_adjust(wspace=0.3)
    
plt.savefig(plotdir + dirname_stripped + '_transport.png', dpi=300)
plt.show()
