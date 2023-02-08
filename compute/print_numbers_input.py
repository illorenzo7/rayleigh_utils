# Author: Loren Matilsky
# Created: 10/21/2022
# print numbers characterizing simulation

import sys, os
sys.path.append(os.environ['raco'])
from cla_util import *
from common import *
from numbers_util import *

# print the purpose of routine
print (buff_line)
print ("NON-DIMENSIONAL NUMBERS (INPUT PARAMETERS)")
print (buff_line)

# read in args
clas0, clas = read_clas(sys.argv)
dirname = clas0.dirname
magnetism = clas0.magnetism
rotation = clas0.rotation

rvals = clas.rvals
if rvals is None:
    rvals = interpret_rvals(dirname, ['rmin', 'rmax'])

for ishell in range(len(rvals) - 1):
    r1 = rvals[ishell]
    r2 = rvals[ishell+1]

    di = get_numbers_input(dirname, r1, r2)
    print (("Shell #%02i: r_1 = " + flt_fmt + ", r_2 = " + flt_fmt)\
            %(ishell + 1, r1, r2))
    print (buff_line)
    lendef1 = 10
    lendef2 = 35
    count = 0
    for key in di.keys():
        print (fill_str(numbers_input_def[key][0], lendef1, ' '), end='')
        print (fill_str(numbers_input_def[key][1], lendef2, ' '), end='')
        print (flt_fmt %di[key])
        if (count + 1) in linebreaks_input:
            print("")
        count += 1
    print (buff_line)
