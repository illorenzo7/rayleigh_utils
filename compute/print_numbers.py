# Author: Loren Matilsky
# Created: 10/21/2022
# print numbers characterizing simulation

import sys, os
sys.path.append(os.environ['raco'])
from cla_util import *
from common import *
from fluid_numbers import *

# read in args
args = sys.argv 
clas0, clas = read_clas(args)
dirname = clas0.dirname
magnetism = clas0.magnetism
rotation = clas0.rotation

rvals = clas.rvals
if rvals is None:
    rvals = interpret_rvals(dirname, ['rmin', 'rmax'])

fmt = "%1.2e"
for ishell in range(len(rvals) - 1):
    r1 = rvals[ishell]
    r2 = rvals[ishell+1]

    di = get_numbers_input(dirname, r1, r2)
    print (buff_line)
    print ("INPUT NON-D NUMBERS")
    print (buff_line)
    print (("Shell #%02i: r1 = " + fmt + ", r2 = " + fmt)\
            %(ishell + 1, r1, r2))
    print (buff_line)
    print (di)
    print (buff_line)
