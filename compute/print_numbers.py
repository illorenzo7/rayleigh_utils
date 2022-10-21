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


