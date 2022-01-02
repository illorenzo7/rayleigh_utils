import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
from common import *
from plotcommon import *
from cla_util import *

# Get CLAs
args = sys.argv 
clas0, clas = read_clas(args)
dirname = clas0.dirname
dirname_stripped = strip_dirname(dirname)

# SPECIFIC ARGS
kwargs_default = dotdict(dict({'the_file': None, 'irvals': np.array([0]), 'rvals': None, 'qvals': np.array([1]), 'modes': [], 'lvals': [], 'mvals': [], 'sectoff': [], 'xminmax': None, 'xmin': None, 'xmax': None, 'yminmax': None, 'ymin': None, 'ymax': None, 'xymin': None, 'xymax': None, 'xyminmax': None}))

# "modes" can be: lpower, mpower, sect, or combinations thereof
# can also get other modes via --lvals (plot freq vs m), --mvals (freq vs l), and --sectoff (m = l - offset)

# many kwargs
find_bad_keys(kwargs_default, clas, 'plot/timetrace/tspec_power_in_range', justwarn=True)

kw = update_dict(kwargs_default, clas)

# get the rvals we want
radlevs = get_slice_levels(dirname)
irvals = kw.irvals
if not kw.rvals is None: # irvals haven't been set directly
    kw.rvals = make_array(kw.rvals)
    if np.all(kw.rvals == 'all'):
        irvals = np.arange(radlevs.nr)
    else:
        irvals = np.zeros_like(kw.rvals, dtype='int')
        for i in range(len(kw.rvals)):
            irvals[i] = np.argmin(np.abs(radlevs.radius/rsun - kw.rvals[i]))

# and the qvals
qvals = make_array(kw.qvals)

# everything must be an array
irvals = make_array(irvals)

# modes, etc.; these must be lists
modes = make_array(kw.modes, tolist=True)
lvals = make_array(kw.lvals, tolist=True)
mvals = make_array(kw.mvals, tolist=True)
sectoff = make_array(kw.sectoff, tolist=True)

everything = modes + lvals + mvals + sectoff
if len(everything) == 0:
    modes = ['lpower']
    everything = modes

print (buff_line)
if kw.the_file is None:
    print ("considering power fractions in temporal spectra")
    print ("qvals = ", qvals)
    print ("irvals = ", irvals)
    nfigures = len(irvals)*len(qvals)*(len(modes) + len(lvals) + len(mvals) + len(sectoff))
else:
    print ("considering ", kw.the_file)
    nfigures = len(modes) + len(lvals) + len(mvals) + len(sectoff)
print ("modes = ", modes)
print ("lvals = ", lvals)
print ("mvals = ", mvals)
print ("sectoff = ", sectoff)
print ("nsets = ", nfigures)
print (buff_line)

# get and make plotdir of non existent
plotdir = my_mkdir(clas0['plotdir'] + 'tspec/')
 
for qval in qvals:
    for irval in irvals:
        # get radial level
        rval = radlevs.radius[irval]/rsun

        # get data
        if kw.the_file is None:
            dataname = ('tspec_qval%04i_irval%02i' %(qval, irval))          
            the_file = get_widest_range_file(clas0['datadir'], dataname)
        else:
            dataname = get_dataname_from_file(kw.the_file)
            the_file = kw.the_file
        iter1, iter2 = get_iters_from_file(the_file)

        #print ("reading " + the_file)
        print ("reading " + the_file)
        di = get_dict(the_file)
        freq = di['freq']
        power_full = np.abs(di['vals'])**2
        nfreq, nell, nm = np.shape(power_full)

        nx, ny = nell, nfreq

        # x and y axes
        x = np.arange(nell)
        y = freq

        # by default sum all the power (should = 1)
        ix1, ix2 = 0, nx - 1
        iy1, iy2 = 0, ny - 1

        # might set both axis boundaries at once with "xy" min/max
        if not kw.xymin is None:
            kw.xmin = kw.xymin
            kw.ymin = kw.xymin
        if not kw.xymax is None:
            kw.xmax = kw.xymax
            kw.ymax = kw.xymax
        if not kw.xyminmax is None:
            kw.xminmax = kw.xyminmax
            kw.yminmax = kw.xyminmax

        # pick and choose part of spectrum to plot
        if not kw.xminmax is None:
            ix1 = np.argmin(np.abs(x - kw.xminmax[0]))
            ix2 = np.argmin(np.abs(x - kw.xminmax[1]))
        if not kw.xmin is None:
            ix1 = np.argmin(np.abs(x - kw.xmin))
        if not kw.xmax is None:
            ix2 = np.argmin(np.abs(x - kw.xmax))

        if not kw.yminmax is None:
            iy1 = np.argmin(np.abs(y - kw.yminmax[0]))
            iy2 = np.argmin(np.abs(y - kw.yminmax[1]))
        if not kw.ymin is None:
            iy1 = np.argmin(np.abs(y - kw.ymin))
        if not kw.ymax is None:
            iy2 = np.argmin(np.abs(y - kw.ymax))

        kw.xminmax = x[ix1], x[ix2] # these are the values we get
        kw.yminmax = y[iy1], y[iy2]

        # add mvals/lvals sectoff to modes
        count = 0
        for mode in everything:
            if count < len(modes):
                basename = mode
                if mode == 'lpower':
                    power = np.sum(power_full, axis=2).T
                if mode == 'mpower':
                    power = np.sum(power_full, axis=1).T
                if mode == 'sect':
                    power = np.zeros((nell, nfreq))
                    for il in range(nell):
                        power[il, :] = power_full[:, il, il]
            elif count < len(modes) + len(lvals):
                lval = int(mode)
                basename = "lval%03i" %lval
                power = power_full[:, lval, :].T
            elif count < len(modes) + len(lvals) + len(mvals):
                mval = int(mode)
                basename = "mval%03i" %mval
                x = np.arange(nell)
                power = power_full[:, :, mval].T
            elif count < len(modes) + len(lvals) + len(mvals) + len(sectoff):
                basename = "sectoff%03i" %mode
                offset = int(mode)
                x = np.arange(nell)
                power = np.zeros((nell, nfreq))
                for il in range(offset, nell):
                    power[il, :] = power_full[:, il, il - offset]

            # now adjust everything by the (x, y) range we want
            power = power[ix1:ix2+1, iy1:iy2+1]

            # Display at terminal what we are analyzing the power of
            if count == 0:
                print ("xminmax = ", kw.xminmax)
                print ("yminmax = ", kw.yminmax)

            print ("mode = " + basename)
            print ("irval = ", irval, '   ', "rval/rsun = %.3f" %rval)

            frac = np.sum(power)/np.sum(power_full)
            print ("frac of total power = " + make_bold("%.3f" %frac))
            print (buff_line)
            count += 1
