# Author: Loren Matilsky
# Created: 10/21/2022
# print output numbers characterizing simulation

import sys, os
sys.path.append(os.environ['raco'])
from cla_util import *
from common import *
from numbers_util import *

# read in args
clas0, clas = read_clas(sys.argv)
dirname = clas0.dirname
dirname_stripped = strip_dirname(dirname)
magnetism = clas0.magnetism
rotation = clas0.rotation

# allowed args + defaults
kwargs_default = dict({'save': False, 'savename': None, 'verbose': False})
kw = update_dict(kwargs_default, clas)

# print the purpose of routine
print (buff_line)
print ('Simulation directory: ', dirname_stripped)
print (buff_line)
print ("NON-DIMENSIONAL NUMBERS (OUTPUT PARAMETERS)")
print (buff_line)

# get desired shells to average over
rvals = clas.rvals
if rvals is None:
    rvals = interpret_rvals(dirname, ['rmin', 'rmax'])

# length of print "tabs" for formatting below
lendef1 = 25
lendef2 = 35

# loop over shells and output numbers
for ishell in range(len(rvals) - 1):
    # print shell info first
    r1 = rvals[ishell]
    r2 = rvals[ishell+1]
    print ('\n\n')
    print (buff_line)
    print (make_bold("Shell #%02i:" %(ishell + 1)))
    print (make_bold('r_1           = ' + flt_fmt %r1))
    print (make_bold('r_2           = ' + flt_fmt %r2))
    print (make_bold('H = r_2 - r_1 = '+ flt_fmt %(r2-r1)))
    print (buff_line)

    # then non-D numbers in shell
    di = get_numbers_output(dirname, r1, r2, verbose=kw.verbose)
    count = 0
    for key in di.keys():
        print (fill_str(numbers_output_def[key][0], lendef1, ' '), end='')
        print (fill_str(numbers_output_def[key][1], lendef2, ' '), end='')
        print (flt_fmt %di[key])
        if (count + 1) in linebreaks_output:
            print("")
        count += 1
    print (buff_line)
