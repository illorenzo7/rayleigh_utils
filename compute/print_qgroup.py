# Author: Loren Matilsky
# Created: 04/26/22
# prints qgroup associated with groupname

import sys, os
sys.path.append(os.environ['raco'])
from cla_util import *
from common import *

# read in args
args = sys.argv
magnetism = False
listall = False

for arg in args:
    if arg == '--mag':
        magnetism = True

groupname = args[1]

try:
    qgroup = get_quantity_group(groupname, magnetism)
except:
    print ('qgroup %s does not exist.' %groupname)
    sys.exit(1)

qvals = qgroup['qvals']

print ('qgroup %s -->' %groupname)
print ("qvals = " + arr_to_str(qvals, "%i"))
