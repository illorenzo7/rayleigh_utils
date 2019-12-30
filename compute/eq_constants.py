# Created: 12/29/2019
# Author: Loren Matilsky

import numpy as np

import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['idref'])

from common import strip_dirname
from reference_tools import equation_coefficients

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

args = sys.argv[2:]
nargs = len(args)
custom_name = None
for i in range(nargs):
    arg = args[i]
    if arg == '-fname':
        custom_name = args[i+1]
    elif arg == '-crb':
        custom_name = 'custom_reference_binary'

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    

eq = equation_coefficients()
if custom_name is None:
    custom_name = 'equation_coefficients'
eq.read(dirname + '/' + custom_name)

print ("For file: ", custom_name)
print ("In directory: ", dirname_stripped)
print ("The equation_coefficients constants are: ")
for i in range(len(eq.constants)): # should always be 10?
    print ('c_%i: %1.3e' %(i + 1, eq.constants[i]))
