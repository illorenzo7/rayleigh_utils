# Created: 11/30/2023
# Author: Loren Matilsky

import numpy as np

import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])

from common import *

jup = False # if True, aspect is considered to be on top (RZ on top), not bottom
args = sys.argv[1:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '--jup':
        jup = True
#    elif arg == '--crb':
#        custom_name = 'custom_reference_binary'

aspect = float(sys.argv[1]) # depth of top layer over bottom layer
base_levels = np.array([0.05, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 0.95])
if jup:
    levsbot = 1./(1.+aspect)*base_levels
    levstop = aspect/(1.+aspect)*base_levels + 1./(1.+aspect)
    levmid = 1./(1.+aspect)
else:
    levsbot = aspect/(1.+aspect)*base_levels
    levstop = 1./(1.+aspect)*base_levels + aspect/(1.+aspect)
    levmid = aspect/(1.+aspect)
all_levels = list(levsbot) + [levmid] + list(levstop)
print(all_levels)
