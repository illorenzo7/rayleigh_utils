# Author: Loren Matilsky
# Created: 10/21/2022
# print input numbers characterizing simulation

import sys, os
sys.path.append(os.environ['raco'])
from cla_util import *
from common import *
from numbers_util import *

# length of print "tabs" for formatting below
lendef1 = 10
lendef2 = 35
buff_line_loc = (lendef1 + lendef2 + 9)*'='

# read in args
clas0, clas = read_clas(sys.argv)
dirname = clas0.dirname
dirname_stripped = strip_dirname(dirname)
magnetism = clas0.magnetism
rotation = clas0.rotation

# allowed args + defaults
kw_default = dict({'save': False, 'savename': None, 'savedir': None, 'verbose': False})
kw = update_dict(kw_default, clas)

# figure out where to potentially save the output
print_funcs = [print]
if kw.save:
    if kw.savename is None:
        kw.savename = 'numbers_input.txt'
    if kw.savedir is None:
        kw.savedir = dirname
    savefile = kw.savedir + '/' + kw.savename
    print (buff_line_loc)
    print ("Saving the output below to:")
    print (savefile)
    f = open(savefile, 'w')
    def my_print(line):
        f.write(line + '\n')
    print_funcs.append(my_print)

# Now print the output, possibly to file as well
for print_func in print_funcs:
    # print the purpose of routine
    print_func (buff_line_loc)
    print_func ('Simulation directory: ' + dirname_stripped)
    print_func ('Full directory: ' + os.path.abspath(dirname))
    print_func (buff_line_loc)
    print_func ("NON-DIMENSIONAL NUMBERS (INPUT PARAMETERS)")
    print_func (buff_line_loc)

    # get desired shells to average over
    rvals = clas.rvals
    if rvals is None:
        rvals = interpret_rvals(dirname, ['rmin', 'rmax'])
    nshells = len(rvals) - 1

    # loop over shells and output numbers
    for ishell in range(nshells):
        # print shell info first
        r1 = rvals[ishell]
        r2 = rvals[ishell+1]
        print_func ('\n\n')
        print_func (buff_line_loc)
        print_func ("Shell %i of %i:" %(ishell + 1, nshells))
        print_func ('r_1           = ' + flt_fmt %r1)
        print_func ('r_2           = ' + flt_fmt %r2)
        print_func ('H = r_2 - r_1 = '+ flt_fmt %(r2-r1))
        print_func (buff_line_loc)

        # then non-D numbers in shell
        di = get_numbers_input(dirname, r1, r2, verbose=kw.verbose)
        count = 0
        for key in di.keys():
            print_func (fill_str(numbers_input_def[key][0], lendef1, ' ') +\
                fill_str(numbers_input_def[key][1], lendef2, ' ') +\
                flt_fmt %di[key])
            if (count + 1) in linebreaks_input:
                print_func("")
            count += 1
        print_func (buff_line_loc)

# remember to close save file
if kw.save:
    f.close()
