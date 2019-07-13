import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

#class ScalarFormatterForceFormat(ScalarFormatter):
#    def _set_format(self,vmin,vmax):  # Override function that finds format to use.
#        self.format = "%1.1f"  # Give format here
#yfmt = ScalarFormatterForceFormat()
#yfmt.set_powerlimits((0,0))

plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import basic_constants as bc

def plotref(r, T, rho, p, dlnT, dlnrho, dlnp, s, dsdr, d2lnrho,\
            label=None, color=None, fig=None, axs=None, ylog=True): 
    figwasNone = False
    if fig is None:     
        figwasNone = True
        fig, axs = plt.subplots(3, 3, figsize= (12,10), sharex=True)
    
    lw = 0.7
    axs[0,0].plot(r/bc.rsun, T, label=label, color=color, linewidth=lw)
#    axs[0,0].yaxis.set_major_formatter(yfmt)
    axs[0,0].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[0,0].set_ylabel(r'$T(r)$' +  ' [K]')
    if ylog:
        axs[0,0].set_yscale('log')
    axs[0,0].set_xlim(np.min(r)/bc.rsun, np.max(r)/bc.rsun)
    
    axs[1,0].plot(r/bc.rsun, rho, label=label, color=color, linewidth=lw)
#    axs[1,0].yaxis.set_major_formatter(yfmt)
    axs[1,0].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[1,0].set_ylabel(r'$\rho(r)\ $' +  r'$\rm{[g\ cm^{-3}]}$')
    if ylog:
        axs[1,0].set_yscale('log')
        
    axs[2,0].plot(r/bc.rsun, p, label=label, color=color, linewidth=lw)
#    axs[2,0].yaxis.set_major_formatter(yfmt)
    axs[2,0].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[2,0].set_xlabel(r'$r/R_\odot$')
    axs[2,0].set_ylabel(r'$p(r)\ $' +  r'$\rm{[erg\ cm^{-3}]}$')
    if ylog:
        axs[2,0].set_yscale('log')
        
    axs[0,1].plot(r/bc.rsun, s, label=label, color=color, linewidth=lw)
#    axs[0,1].yaxis.set_major_formatter(yfmt)
    axs[0,1].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[0,1].set_xlabel(r'$r/R_\odot$')
    axs[0,1].set_ylabel(r'$s(r)\ $' +  r'$\rm{[erg\ g^{-1}\ K^{-1}]}$') 
        
    axs[1,1].plot(r/bc.rsun, dlnT, label=label, color=color, linewidth=lw)
#    axs[1,1].yaxis.set_major_formatter(yfmt)
    axs[1,1].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[1,1].set_ylabel(r'$d\ln{T}/dr\ $' +  r'$\rm{[cm^{-1}]}$')
        
    axs[2,1].plot(r/bc.rsun, dlnrho, label=label, color=color, linewidth=lw)
#    axs[2,1].yaxis.set_major_formatter(yfmt)
    axs[2,1].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[2,1].set_ylabel(r'$d\ln{\rho}/dr\ $' +  r'$\rm{[cm^{-1}]}$')
        
    axs[0,2].plot(r/bc.rsun, dlnp, label=label, color=color, linewidth=lw)
#    axs[0,2].yaxis.set_major_formatter(yfmt)
    axs[0,2].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[0,2].set_xlabel(r'$r/R_\odot$')
    axs[0,2].set_ylabel(r'$d\ln{p}/dr\ $' +  r'$\rm{[cm^{-1}]}$')  

    axs[1,2].plot(r/bc.rsun, dsdr, label=label, color=color, linewidth=lw)
#    axs[1,2].yaxis.set_major_formatter(yfmt)
    axs[1,2].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[1,2].set_xlabel(r'$r/R_\odot$')
    axs[1,2].set_ylabel(r'$ds/dr\ $' +\
       r'$\rm{[erg\ g^{-1}\ K^{-1}\ cm^{-1}]}$')  
        
    axs[2,2].plot(r/bc.rsun, d2lnrho, label=label, color=color, linewidth=lw)
#    axs[2,2].yaxis.set_major_formatter(yfmt)
    axs[2,2].ticklabel_format(scilimits = (0,0), useMathText=True, axis='y')
    axs[2,2].set_xlabel(r'$r/R_\odot$')
    axs[2,2].set_ylabel(r'$d^2\ln{\rho}/dr^2\ $' +\
       r'$\rm{[cm^{-2}]}$')  
        
    # Get ticks everywhere
    plt.sca(axs[0,0])
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')
    
    plt.sca(axs[1,0])
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')
    
    plt.sca(axs[2,0])
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')
 
    plt.sca(axs[0,1])
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')    
    
    plt.sca(axs[1,1])
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')  
    
    plt.sca(axs[2,1])
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')
    
    plt.sca(axs[0,2])
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both') 
    
    plt.sca(axs[1,2])
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')
    
    plt.sca(axs[2,2])
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')
    
    if figwasNone:
        return fig, axs
