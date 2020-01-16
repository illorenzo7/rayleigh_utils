###############################################
# Author: Loren Matilsky
# Date created: 01/13/2020
#
# This script plots takes radial energy fluxes from Shell_Avgs,
# and the existing heating function, and modifies the heating function
# to bring the system into flux balance

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.interpolate import interp1d
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import get_widest_range_file, strip_dirname,\
        get_iters_from_file, get_dict, rsun
from get_parameter import get_parameter
from reference_tools import equation_coefficients

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Find the Shell_Avgs file(s) in the data directory. If there are multiple, 
# by default choose the one with widest range in the average
Shell_Avgs_file = get_widest_range_file(datadir, 'Shell_Avgs')

# Get the original reference state with fine grid
eq = equation_coefficients()
eq.read(dirname + '/custom_reference_binary')
nfine = eq.nr
rr_fine = eq.radius
print ("custom_reference_binary nr: %i" %nfine)

# Get command-line arguments to adjust the interval of averaging files
lum = 3.846e33
plot_chop = False
window_width = int(0.02*nfine)

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-usefile':
        Shell_Avgs_file = args[i+1]
        Shell_Avgs_file = Shell_Avgs_file.split('/')[-1]
    elif arg == '-lum':
        lum = float(args[i+1])
    elif arg == '-chop':
        plot_chop = True
    elif arg == '-window':
        window_width = int(args[i+1])

#Create the plot
lw = 1. # regular lines
#lw = 1.5 # Bit thicker lines

# Read in the flux data
print ('Getting radial fluxes from ' + datadir + Shell_Avgs_file + ' ...')
di = get_dict(datadir + Shell_Avgs_file)
vals = di['vals']
lut = di['lut']
nq = di['nq']
rr = di['rr']
nr = di['nr']
ri = di['ri']
ro = di['ro']

# Determine the simulation is magnetic
magnetism = get_parameter(dirname, 'magnetism')

qindex_hflux = lut[1433]
qindex_cflux = lut[1470]
qindex_kflux = lut[1923]
qindex_vflux = lut[1935]
qindex_eflux = lut[1455]

hflux = vals[:, lut[1433]]
eflux = vals[:, lut[1455]]
cflux = vals[:, lut[1470]]
kflux = vals[:, lut[1923]]
vflux = -vals[:, lut[1935]]
tflux = hflux + eflux + cflux + kflux + vflux 

if magnetism:
    qindex_mflux = lut[2001] # this is actually (-4*pi) TIMES 
                        # the correct Poynting flux
    mflux = -vals[:, qindex_mflux]/(4*np.pi)
    tflux += mflux


# Interpolate the total flux onto a fine grid
func = interp1d(rr, tflux, fill_value='extrapolate')
tflux_fine = func(rr_fine)

# Calculate the heat flux needed to add
tflux_eq = lum/(4*np.pi*rr_fine**2)
print ("Computing change to heating function for lum = %1.3e" %lum)
q = 1/rr_fine**2*np.gradient(rr_fine**2*tflux_fine, rr_fine)

# Pad q with zeros at the ends (don't want to change the heating here anyway
# )
nzero = int(0.05*nfine)
print("Padding q(r) with %i zeros at each end" %nzero)
q[:nzero] = 0.0
q[-nzero:] = 0.0

# Smooth q box window function (sliding mean)
print ("Smoothing q with window_width = %i" %window_width)
print ("To specify alternate width, use -window option")
q_smooth = np.zeros(nfine)
for ir in range(nfine):
    ir_min = max(0, ir - window_width//2)
    ir_max = min(nfine - 1, ir + window_width//2)
    q_smooth[ir] = np.mean(q[ir_min:ir_max + 1])

fig,axs = plt.subplots(1,2)
plt.sca(axs[0])
plt.plot(rr, tflux, 'b', label='F_tot, orig')
plt.plot(rr_fine, tflux_fine, 'r--', label='F_tot, fine interp')
plt.plot(rr_fine, tflux_eq, 'g', label='F_tot, eq')
plt.legend()
plt.sca(axs[1])
if plot_chop:
    plt.plot(rr_fine, q, 'b', label='choppy q')
plt.plot(rr_fine, q_smooth, 'r--', label='smooth q')
plt.legend()
plt.show()

Q_orig = eq.constants[9]*eq.functions[5]
Q_new = Q_orig + q_smooth

ir_min = np.argmin(np.abs(rr_fine - ri))
ir_max = np.argmin(np.abs(rr_fine - ro))
int_profile = -4*np.pi*simps(Q_new[ir_max:ir_min + 1]*\
        rr_fine[ir_max:ir_min + 1]**2.0, rr_fine[ir_max:ir_min + 1]) # remember rr is reversed
print("Integral of Q(r) + q(r) is %1.3e" %int_profile)
print("Now renormalizing Q(r) + q(r)")
radial_shape = Q_new/int_profile

eq.set_function(radial_shape, 6)
savename = 'custom_reference_binary2'
print ("Saving new reference state (with mod heating) to ", savename)
eq.write(dirname + '/' + savename)
