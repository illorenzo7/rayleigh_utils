# Author: Loren Matilsky
# Created: 09/11/2020
#
# Purpose: generate a magnetic reference file from a hydro reference file

# Parameters: dirname (first argument), 
# Command-line options:
#
import numpy as np
import sys, os, shutil

sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])

from reference_tools import equation_coefficients
from common import *
import basic_constants as bc

# Get directory name
dirname = sys.argv[1]
fname = 'custom_reference_binary'

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-fname':
        polytropic_reference = False
        fname = args[i+1]

# Back up current reference file and read in reference state
print('backing up reference file in ' + dirname + '/:')
fname_bak = fname + '_bak'
print('%s --> %s' %(fname, fname_bak))
shutil.copyfile(dirname + '/' + fname, dirname + '/' + fname_bak)
print('Reading reference file ' + dirname + '/' + fname)
eq = equation_coefficients()
eq.read(dirname + '/' + fname)

# Set the magnetic constants in the original file
print('setting c_4 = c_9 = 1/(4*pi), c_7  = 1 in %s' %fname)
eq.set_constant(1.0/4.0/np.pi, 4) # multiplies Lorentz force
eq.set_constant(1.0, 7) # multiplies eta in induction-diffusion term
eq.set_constant(1.0/4.0/np.pi, 9) # multiplies Joule heating
print('writing %s' %fname)
eq.write(dirname + '/' + fname)
