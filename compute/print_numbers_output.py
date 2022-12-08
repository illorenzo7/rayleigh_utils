# Author: Loren Matilsky
# Created: 10/21/2022
# print numbers characterizing simulation

import sys, os
sys.path.append(os.environ['raco'])
from cla_util import *
from common import *
from fluid_numbers import *

# print the purpose of routine
print (buff_line)
print ("NON-DIMENSIONAL NUMBERS (OUTPUT PARAMETERS)")
print (buff_line)

# read in args
args = sys.argv 
clas0, clas = read_clas(args)
dirname = clas0.dirname
magnetism = clas0.magnetism
rotation = clas0.rotation

# allowed args + defaults, then update
kwargs_default = dict({'the_file': None, 'the_file_az': None, 'sd': None})
kw = update_dict(kwargs_default, clas)

# need these likely
rmin, rmax = interpret_rvals(dirname, ['rmin', 'rmax'])

# deal with shell depth (by default use whole shell)
shell_depth = clas.sd
if shell_depth is None:
    shell_depth = rmax - rmin

# print the shell depth we use
fmt = "%1.2e"
print (("System shell depth: H = " + fmt + " cm") %shell_depth)

# deal with rvals to average numbers over
rvals = clas.rvals
if rvals is None:
    rvals = rmin, rmax

# get the output numbers
di = get_numbers_output(dirname, shell_depth, kw.the_file, kw.the_file_az)

# now print the numbers
for ishell in range(len(rvals) - 1):
    r1 = rvals[ishell]
    r2 = rvals[ishell+1]

    print (("Shell #%02i: r_1 = " + fmt + ", r_2 = " + fmt)\
            %(ishell + 1, r1, r2))
    print (buff_line)
    lendef1 = 25
    lendef2 = 35
    count = 0
    for key in di.keys():
        print (fill_str(numbers_output_def[key][0], lendef1, ' '), end='')
        print (fill_str(numbers_output_def[key][1], lendef2, ' '), end='')
        # get the actual number (volume average)
        num = volav_in_radius(dirname, di[key], r1, r2)
        print (fmt %num)
        if (count + 1) in linebreaks_output:
            print("")
        count += 1
    print (buff_line)
