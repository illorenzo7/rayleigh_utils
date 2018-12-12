# Author Loren Matilsky
# Date created: 06/08/2017

#  This example routine makes use of the GlobalAverage
#  data structure associated with the G_Avg output.
#  Upon initializing a GlobalAverage object, the 
#  object will contain the following attributes:
#
#    ----------------------------------
#    self.nrec                  : number of time steps
#    self.nq                     : number of diagnostic quantities output
#    self.qv[0:nq-1]             : quantity codes for the diagnostics output
#    self.vals[0:nrec-1,0:nq-1] : Globally averaged diagnostics as function of time and quantity index
#    self.iters[0:nrec-1]       : The time step numbers stored in this output file
#    self.time[0:nrec-1]        : The simulation time corresponding to each time step
#    self.lut                    : Lookup table for the different diagnostics output

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys, os
from subprocess import call
from common import get_file_lists, get_widest_range_file

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]

# Data and plot directories
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)
    
# Find the etrace file(s) in the data directory. If there are multiple, by
# default choose the one with widest range in the trace.
etrace_file = get_widest_range_file(datadir, 'etrace')

# Set defaults
showplot = False
xiter = False
notfrom0 = False

# Get command-line arguments
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-show'):
        showplot = True
    elif (arg == '-xiter'): # plot w.r.t. iterations
        xiter = True
    elif (arg == '-usefile'):
        etrace_file = args[i+1]
        etrace_file = etrace_file.split('/')[-1]
    elif (arg == '-notfrom0'):
        notfrom0 = True

# Get the stripped etrace file to use in the plot name
etrace_file_stripped = etrace_file.split('.')[0]
if (xiter):
    tag = '_xiter'
else:
    tag = '_xtime'
savename = etrace_file_stripped + tag + '.png'

# Read in the KE data
print ('Reading in KE components data from ' + datadir + etrace_file + ' ...')
aaaa

if (not xiter):
    xaxis = times
else:
    xaxis = iters

if (notfrom0):
    x_min = np.min(xaxis)
else:
    x_min = 0

# create figure
fig=plt.figure(1,figsize=(15,10))

# first plot: total kinetic energy trace
ax1 = fig.add_subplot(221)
# actual plotting
ax1.plot(xaxis,ke,'black', label = r'$\rm{KE_{total}}$'+\
        r'$ \  = \ \frac{1}{2}\overline{\rho} v^2$')
ax1.plot(xaxis,rke,'r', label = r'$\rm{KE}$' + \
        r'$_r \  = \ \frac{1}{2}\overline{\rho} v_r^2$')
ax1.plot(xaxis,tke,'g', label = r'$\rm{KE}$'+ \
        r'$_\theta \  = \ \frac{1}{2}\overline{\rho} v_\theta^2$')
ax1.plot(xaxis,pke,'b', label = r'$\rm{KE}$'+\
        r'$_\phi \  = \ \frac{1}{2}\overline{\rho} v_\phi^2$')
# yscale
#ax1.set_yscale('log')
# title and axis labels
if (xiter):
    ax1.set_xlabel('iteration #',fontsize=14)
else:
    ax1.set_xlabel('t (days)',fontsize=14)
ax1.set_ylabel(r'$\rm{Energy} \ \rm{Density}\ ( \rm{erg}\ \rm{cm}^{-3})$',\
        fontsize=14)
ax1.set_title('Global average of total KE', fontsize=16)
# Set limits
ax1.set_xlim((x_min,np.max(xaxis)))
ax1.ticklabel_format(scilimits = (0,0), useMathText=True)
#minval1 = min(np.min(rke),np.min(tke),np.min(pke))
maxval1 = np.max(ke) 
ax1.set_ylim((0.,maxval1*1.05))
# legend
#ax1.legend(bbox_to_anchor=(1.2,0.75),shadow=True) 
ax1.legend()

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')

# Make the second plot (kinetic energy of the mean motions)
ax2 = fig.add_subplot(222)
# actual plotting
ax2.plot(xaxis,mke,'k', label =\
        r'$\frac{1}{2}\overline{\rho}\langle |\mathbf{v}|\rangle^2$')
ax2.plot(xaxis,mrke,'r', \
        label = r'$\frac{1}{2}\overline{\rho}\langle v_r\rangle^2$')
ax2.plot(xaxis,mtke,'g',\
        label = r'$\frac{1}{2}\overline{\rho}\langle v_\theta\rangle^2$')
ax2.plot(xaxis,mpke,'b',\
        label = r'$\frac{1}{2}\overline{\rho}\langle v_\phi\rangle^2$')
# yscale
ax2.set_yscale('log')
# set limits
ax2.set_xlim((x_min, np.max(xaxis)))
minval2 = min(np.min(mrke), np.min(mtke), np.min(mpke))
maxval2 = np.max(mke)
ax2.set_ylim((minval2/3., maxval2*3.))
# title and axis labels
ax2.set_title('mean KE', fontsize=16)
if (xiter):
    ax2.set_xlabel('iteration #',fontsize=14)
else:
    ax2.set_xlabel('t (days)',fontsize=14)
ax2.set_ylabel(r'$\rm{Energy} \ \rm{Density}\ ( \rm{erg}\ \rm{cm}^{-3})$',\
        fontsize=14)
ax2.ticklabel_format(axis='x', scilimits = (0,0), useMathText=True)
# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')

#ax2.ticklabel_format(scilimits = (0,0), useMathText=True)

# Make the third plot (kinetic energy of the fluctuating, or convective,
        # motions
ax3 = fig.add_subplot(223)
# actual plotting
ax3.plot(xaxis,fke,'k', label =\
        r'$\frac{1}{2}\overline{\rho}\langle |\mathbf{v}^\prime|^2\rangle$')
ax3.plot(xaxis,frke,'r', label =\
        r'$\frac{1}{2}\overline{\rho}\langle (v_r^\prime)^2\rangle$')
ax3.plot(xaxis,ftke,'g', label =\
        r'$\frac{1}{2}\overline{\rho}\langle (v_\theta^\prime)^2\rangle$')
ax3.plot(xaxis,fpke,'b', label =\
        r'$\frac{1}{2}\overline{\rho}\langle v_\phi^\prime)^2\rangle$')
# set limits
ax3.set_xlim((x_min, np.max(xaxis)))
minval3 = min(np.min(frke), np.min(ftke), np.min(mpke))
maxval3 = np.max(fke)
ax3.set_ylim((0., maxval3*1.05))
# title and axis labels
ax3.set_title(' fluc KE', fontsize=16)
if (xiter):
    ax3.set_xlabel('iteration #',fontsize=14)
else:
    ax3.set_xlabel('t (days)',fontsize=14)
ax3.set_ylabel(r'$\rm{Energy} \ \rm{Density}\ ( \rm{erg}\ \rm{cm}^{-3})$',\
        fontsize=14)
ax3.ticklabel_format(scilimits = (0,0), useMathText=True)

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')

# Make the fourth plot (sum of mean and fluctuating KE, just to check that
# it equals the total
tke = fke + mke
trke = frke + mrke
ttke = ftke + mtke
tpke = fpke + mpke

ax4 = fig.add_subplot(224)
# actual plotting
ax4.plot(xaxis,tke,'k', label = 'mean + fluc (t)')
ax4.plot(xaxis,trke,'r', label = 'mean + fluc (r)')
ax4.plot(xaxis,ttke,'g', label = 'mean + fluc (t)')
ax4.plot(xaxis,tpke,'b', label = 'mean + fluc (p)')
# set limits
ax4.set_xlim((x_min, np.max(xaxis)))
maxval4 = np.max(tke)
ax4.set_ylim((0., maxval4*1.05))
# title and axis labels
ax4.set_title(' mean + fluc KE', fontsize=16)
if (xiter):
    ax4.set_xlabel('iteration #',fontsize=14)
else:
    ax4.set_xlabel('t (days)',fontsize=14)
ax4.set_ylabel(r'$\rm{Energy} \ \rm{Density}\ ( \rm{erg}\ \rm{cm}^{-3})$',\
        fontsize=14)
ax4.ticklabel_format(scilimits = (0,0), useMathText=True)

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')

# Space the subplots to make them look pretty
plt.tight_layout
#plt.subplots_adjust(right=0.92, left=0.08, bottom=0.1, top=0.90,
#        wspace=0.4, hspace=0.4)
plt.subplots_adjust(right=0.92, left=0.08, bottom=0.1, top=0.90,
        hspace=0.4)

# Save the plot
print ('Saving the etrace plot at ' + plotdir + savename + ' ...')
plt.savefig(plotdir + savename)
plt.close()

# Call 'eog' to show the plot if demanded by the user
if (showplot):
    call (['eog', plotdir + savename])
