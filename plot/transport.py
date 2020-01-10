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
from rayleigh_diagnostics import TransportCoeffs
from reference_tools import equation_coefficients
from read_kappa00 import read_kappa00

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
fname = None
xminmax = None
plot_kappa00 = False
for i in range(nargs):
    arg = args[i]
    if arg == '-log':
        ylog = True
    elif arg == '-fname':
        fname = args[i+1]
    elif arg == '-crb':
        fname = 'custom_reference_binary'
    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-00':
        plot_kappa00 = True

magnetism = get_parameter(dirname, 'magnetism')

if not fname is None:
    print ("Getting transport coefs from ", fname)
    eq = equation_coefficients()
    eq.read(dirname + '/' + fname)
    r = eq.radius
    nu = eq.constants[4]*eq.functions[2, :]
    dlnu = eq.functions[10, :]
    kappa = eq.constants[5]*eq.functions[4, :]
    dlnkappa = eq.functions[11, :]
    if magnetism:
        eta = eq.constants[6]*eq.functions[6, :]
        dlneta = eq.functions[12, :]
else:
    try:
        t = TransportCoeffs(dirname + '/transport', '')
        print ("Getting transport coefs from 'transport' file")
        r = t.radius
        nu = t.nu
        dlnu = t.dlnu
        kappa = t.kappa
        dlnkappa = t.dlnkappa
        if magnetism:
            eta = t.eta
            dlneta = t.dlneta
    except:
        print ("Getting transport coefs from 'equation_coefficients' file")
        eq = equation_coefficients()
        eq.read(dirname + '/equation_coefficients')
        r = eq.radius
        nu = eq.constants[4]*eq.functions[2, :]
        dlnu = eq.functions[10, :]
        kappa = eq.constants[5]*eq.functions[4, :]
        dlnkappa = eq.functions[11, :]
        if magnetism:
            eta = eq.constants[6]*eq.functions[6, :]
            dlneta = eq.functions[12, :]

if plot_kappa00:
    print ("Getting kappa00 coefficient from 'kappa00' file")
    dummy, dummy, kappa00, dlnkappa00 = read_kappa00(dirname + '/kappa00') 
    # nr and rr are "dummies"

if magnetism:
    if plot_kappa00:
        figsize = 7., 13.
        fig, axs = plt.subplots(4, 2, figsize=figsize, sharex=True)
    else:
        figsize = 7., 10.
        fig, axs = plt.subplots(3, 2, figsize=figsize, sharex=True)
else:
    if plot_kappa00:
        figsize = 7., 10.
        fig, axs = plt.subplots(3, 2, figsize=figsize, sharex=True)
    else:
        figsize = 7., 7.
        fig, axs = plt.subplots(2, 2, figsize=figsize, sharex=True)

# If xminmax was specified only a plot a portion of the profiles
if not xminmax is None:
    if r[0] < r[-1]: # radii in increasing order
        ir1 = np.argmin(np.abs(r - xminmax[0]))
        ir2 = np.argmin(np.abs(r - xminmax[1]))
    else: # radii in decreasing order
        ir1 = np.argmin(np.abs(r - xminmax[1]))
        ir2 = np.argmin(np.abs(r - xminmax[0]))
    r = r[ir1:ir2+1]
    nu = nu[ir1:ir2+1]
    kappa = kappa[ir1:ir2+1]
    dlnu = dlnu[ir1:ir2+1]
    dlnkappa = dlnkappa[ir1:ir2+1]
    if magnetism:
        eta = eta[ir1:ir2+1]
        dlneta = dlneta[ir1:ir2+1]
    if plot_kappa00:
        kappa00 = kappa00[ir1:ir2+1]
        dlnkappa00 = dlnkappa00[ir1:ir2+1]

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

mag_index = 2
if plot_kappa00:
    mag_index += 1
    axs[2,0].plot(r/rsun, kappa00, linewidth=lw)
    #    axs[1,0].yaxis.set_major_formatter(yfmt)
    axs[2,0].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[2,0].set_ylabel(r'$\kappa_{00}(r)\ $' +  r'$\rm{[cm^2\ s^{-1}]}$')
    if ylog:
        axs[2,0].set_yscale('log')

    axs[2,1].plot(r/rsun, dlnkappa00, linewidth=lw)
    #    axs[1,0].yaxis.set_major_formatter(yfmt)
    axs[2,1].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[2,1].set_ylabel(r'$d\ln\kappa_{00}/dr\ $' +  r'$\rm{[cm^{-1}]}$')
    
    # ticks everywhere
    plt.sca(axs[2,0])
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    plt.sca(axs[2,1])
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')    

if magnetism:
    axs[mag_index,0].plot(r/rsun, eta, linewidth=lw)
    #    axs[1,0].yaxis.set_major_formatter(yfmt)
    axs[mag_index,0].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[mag_index,0].set_ylabel(r'$\eta(r)\ $' +  r'$\rm{[cm^2\ s^{-1}]}$')
    if ylog:
        axs[mag_index,0].set_yscale('log')

    axs[mag_index,1].plot(r/rsun, dlneta, linewidth=lw)
    #    axs[1,0].yaxis.set_major_formatter(yfmt)
    axs[mag_index,1].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[mag_index,1].set_ylabel(r'$d\ln\eta/dr\ $' +  r'$\rm{[cm^{-1}]}$')
    
    # ticks everywhere
    plt.sca(axs[mag_index,0])
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    plt.sca(axs[mag_index,1])
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
