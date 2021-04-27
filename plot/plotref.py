import sys, os
import numpy as np
import matplotlib.pyplot as plt
sys.path.append(os.environ['raco'] + '/tachocline')
from common import *

plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}

def plotref(r, T, rho, p, dlnT, dlnrho, s, dsdr, d2lnrho,\
        gravity, Q,\
            label=None, color=None, fig=None, axs=None, ylog=True,\
            xminmax=None): 
    if not xminmax is None:
        ir_min = np.argmin(np.abs(r - xminmax[0]))
        ir_max = np.argmin(np.abs(r - xminmax[1]))
        r = np.copy(r[ir_max:ir_min + 1])
        T = np.copy(T[ir_max:ir_min + 1])
        rho = np.copy(rho[ir_max:ir_min + 1])
        p = np.copy(p[ir_max:ir_min + 1])
        dlnT = np.copy(dlnT[ir_max:ir_min + 1])
        dlnrho = np.copy(dlnrho[ir_max:ir_min + 1])
        s = np.copy(s[ir_max:ir_min + 1])
        dsdr = np.copy(dsdr[ir_max:ir_min + 1])
        d2lnrho = np.copy(d2lnrho[ir_max:ir_min + 1])
        gravity = np.copy(gravity[ir_max:ir_min + 1])
        Q = np.copy(Q[ir_max:ir_min + 1])

    figwasNone = False
    if fig is None:     
        figwasNone = True
        fig, axs = plt.subplots(3, 4, figsize= (10,6), sharex=True)
    
    lw = 0.7
    axs[0,0].plot(r/rsun, T, label=label, color=color, linewidth=lw)
    axs[0,0].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[0,0].set_ylabel(r'$T(r)$' +  ' [K]')
    if ylog:
        axs[0,0].set_yscale('log')
    axs[0,0].set_xlim(np.min(r)/rsun, np.max(r)/rsun)
    
    axs[0,1].plot(r/rsun, rho, label=label, color=color, linewidth=lw)
    axs[0,1].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[0,1].set_ylabel(r'$\rho(r)\ $' +  r'$\rm{[g\ cm^{-3}]}$')
    if ylog:
        axs[0,1].set_yscale('log')
        
    axs[0,2].plot(r/rsun, p, label=label, color=color, linewidth=lw)
    axs[0,2].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[0,2].set_xlabel(r'$r/R_\odot$')
    axs[0,2].set_ylabel(r'$p(r)\ $' +  r'$\rm{[erg\ cm^{-3}]}$')
    if ylog:
        axs[0,2].set_yscale('log')
        
    axs[0,3].plot(r/rsun, s, label=label, color=color, linewidth=lw)
    axs[0,3].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[0,3].set_xlabel(r'$r/R_\odot$')
    axs[0,3].set_ylabel(r'$s(r)\ $' +  r'$\rm{[erg\ g^{-1}\ K^{-1}]}$') 
        
    axs[1,0].plot(r/rsun, dlnT, label=label, color=color, linewidth=lw)
    axs[1,0].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[1,0].set_ylabel(r'$d\ln{T}/dr\ $' +  r'$\rm{[cm^{-1}]}$')
        
    axs[1,1].plot(r/rsun, dlnrho, label=label, color=color, linewidth=lw)
    axs[1,1].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[1,1].set_ylabel(r'$d\ln{\rho}/dr\ $' +  r'$\rm{[cm^{-1}]}$')
        
    axs[1,2].plot(r/rsun, d2lnrho, label=label, color=color, linewidth=lw)
    axs[1,2].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[1,2].set_xlabel(r'$r/R_\odot$')
    axs[1,2].set_ylabel(r'$d^2\ln{\rho}/dr^2\ $' +\
       r'$\rm{[cm^{-2}]}$')  

    axs[1,3].plot(r/rsun, dsdr, label=label, color=color, linewidth=lw)
    axs[1,3].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[1,3].set_xlabel(r'$r/R_\odot$')
    axs[1,3].set_ylabel(r'$ds/dr\ $' +\
       r'$\rm{[erg\ g^{-1}\ K^{-1}\ cm^{-1}]}$')  

    # Gravity and heating, last row
    axs[2,0].plot(r/rsun, gravity, label=label, color=color, linewidth=lw)
    axs[2,0].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[2,0].set_xlabel(r'$r/R_\odot$')
    axs[2,0].set_ylabel(r'$g(r)\ $' +\
       r'$\rm{[cm\ s^{-2}$')  

    axs[2,1].plot(r/rsun, Q, label=label, color=color, linewidth=lw)
    axs[2,1].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[2,1].set_xlabel(r'$r/R_\odot$')
    axs[2,1].set_ylabel(r'$Q(r)\ $' + r'$\rm{[erg\ cm^{-3}\ s^{-1}]}$')  
      
    # Get ticks everywhere
    for ax in axs.flatten():
        plt.sca(ax)
        plt.minorticks_on()
        plt.tick_params(top=True, right=True, direction='in', which='both')
    
    if figwasNone:
        return fig, axs
