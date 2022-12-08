import numpy as np
import matplotlib.pyplot as plt
import sys
import os

dirname = sys.argv[1]
radatadir = dirname + '/Shell_Avgs/'

datadir = dirname + '/data/'
#plotdir = dirname + '/plots/'
plotdir = '/altair/loma3853/rayleigh/n3_plots/s_az/'
dirname_stripped = dirname.split('/')[-1]

s_sph = np.load(datadir + 's_spherical_mean.npy')
s_az = np.load(datadir + 's_azimuthal_mean.npy')

s_fluc = s_az - s_sph.reshape((1,128))