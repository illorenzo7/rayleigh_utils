# Routine to trace Rayleigh G_Avgs data in time
# Created by: Loren Matilsky
# On: 03/11/2020
############################################################################
# This routine computes a trace in time/spherical harmonics l and m
# contained in file [fname] (produced by compute/timetrace/lmvals_v) and
# computes the power spectrum by doing a Fourier
# transform in time
# spectrum is 3D (function of l, m and omega)
# 
# stores data in a new .pkl file, named after the lmvals_v trace:
# [dirname]_trace_lmvals_om_[iter1]_[iter2].pkl -->
# [dirname]_lmvals_om_pspec_[iter1]_[iter2].pkl
#
# The final datacube output ('vals') will have shape
# (nfreq, nell, ndepths, 3) = (ntimes, nell, ndepths, 3)

# Import relevant modules
import numpy as np
import pickle
import sys, os
sys.path.append(os.environ['raco'])
from common import *

# File name of the trace_lmvals_v file 
fname = sys.argv[1]

# New file replacing trace_lmvals_om --> lmvals_om_pspec
try:
    savename = fname.replace('trace_lmvals_om', 'lmvals_om_pspec')
except:
    print("ERROR: input file must be of type 'trace_lmvals_om'")
    print("Exiting")
    sys.exit()

# Read in trace_lmvals_om
print("Reading ", fname)
di_in = get_dict(fname)

vals_omr = di_in['vals_omr']
vals_omt = di_in['vals_omt']
vals_omp = di_in['vals_omp']
times = di_in['times']
iters = di_in['iters']
iter1 = di_in['iter1']
iter2 = di_in['iter2']
niter = di_in['niter']
lut = di_in['lut']
qv = di_in['qv']
nq = di_in['nq']
rvals = di_in['rvals']
rinds = di_in['rinds']
nr = di_in['nr']
lvals = di_in['lvals']
mvals = di_in['mvals']
nell = di_in['nell']
lmax = di_in['lmax']

# Check to see if time-samples are equally spaced
dtimes = np.diff(times)
dt = np.mean(dtimes)
dt_spread = np.std(dtimes)
if dt_spread == 0.:
    print("check: samples are equally spaced in time")
else:
    print("WARNING! samples not equally spaced")
    print("sigma(dtimes)/mean(dtimes) = ", dt_spread/dt)

# Compute the Nyquist frequency and frequency spacing
f_nyq = 1./2./dt
df = 1./niter/dt

# Do the FFT along time-axis
print("Performing FFT")
omr_fft = np.abs(np.fft.fft(vals_omr, axis=0))**2
omt_fft = np.abs(np.fft.fft(vals_omt, axis=0))**2
omp_fft = np.abs(np.fft.fft(vals_omp, axis=0))**2
freq = np.fft.fftfreq(niter, dt)

print("Shifting the array so zero-frequency is in center")
omr_fft = np.fft.fftshift(omr_fft, axes=0)
omt_fft = np.fft.fftshift(omt_fft, axes=0)
omp_fft = np.fft.fftshift(omp_fft, axes=0)
freq = np.fft.fftshift(freq)
nfreq = len(freq)

# ND-array of frequencies
freq_nd = freq.reshape((nfreq, 1, 1, 1))

# Save the power spectrum
print ('Saving power spectrum at ' + savename)
f = open(savename, 'wb')
pickle.dump({'omr_fft': omr_fft, 'omt_fft': omt_fft, 'omp_fft': omp_fft,\
        'freq': freq, 'freq_nd': freq_nd, 'f_nyq': f_nyq, 'df': df,\
        'dt': dt, 'dt_spread': dt_spread,\
        'times': times, 'iters': iters, 'iter1': iter1, 'iter2': iter2,\
        'niter': niter, 'lut': lut, 'qv': qv, 'nq': nq, 'rvals': rvals,\
        'rinds': rinds, 'nr': nr, 'lvals': lvals, 'mvals': mvals,\
        'nell': nell, 'lmax': lmax}, f, protocol=4)
f.close()
