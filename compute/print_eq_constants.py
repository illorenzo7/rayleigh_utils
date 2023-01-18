# Created: 12/29/2019
# Author: Loren Matilsky

import numpy as np

import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])

from common import *
from reference_tools import equation_coefficients

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

args = sys.argv[2:]
nargs = len(args)
custom_name = None
save = False
for i in range(nargs):
    arg = args[i]
    if arg == '--fname':
        custom_name = args[i+1]
    elif arg == '--crb':
        custom_name = 'custom_reference_binary'
    elif arg == '--save':
        save = True

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    

eq = equation_coefficients()
c_dict = reverse_dict(eq.c_dict)

if custom_name is None:
    custom_name = 'equation_coefficients'
eq.read(dirname + '/' + custom_name)

print (buff_line)
print ("For file: ", custom_name)
print ("In directory: ", dirname_stripped)
print ("The equation_coefficients constants are: ")
if save:
    savename = dirname + '/eq_constants.txt'
    f = open(savename, 'w')
for i in range(len(eq.constants)): # should always be 10?
    cnum = fill_str('c_%i' %(i+1), 4, ' ')
    cname = fill_str(c_dict[i+1], 15, ' ')
    print ('%s = %s = %1.3e' %(cnum, cname, eq.constants[i]))
    if save:
        f.write('%s = %s = %1.7e\n' %(cnum, cname, eq.constants[i]))
if save:
    f.close()
    print (buff_line)
    print ("saving constants in", savename)
    print (buff_line)
