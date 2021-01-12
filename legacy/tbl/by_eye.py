# Author: Loren Matilsky
# Created on: 05/04/2018
# Computes the radius at the base of the thermal boundary layer using input 
# from user, who determines this depth "bye eye" (most likely by examining
# the conductive flux profile with radius)
# Usage: python compute/by_eye.py [dirname] [r_tbl] [[-ro]]
# saves radial index of r_tbl  at data/[dirname]_ir_tbl.npy

# Import relevant modules here
import numpy as np
import sys

# Get simulation directory and data directory
dirname = sys.argv[1]
dirname_stripped = dirname.split('/')[-1]
datadir = dirname + '/data/'

# Get desired location of the TBL
r_tbl = float(sys.argv[2])

args = sys.argv[3:]
nargs = len(args)

# by default, treat radius user enters as modulated by the solar radius
normalize_by_ro = False
for i in range(nargs):
    arg = args[i]
    if (arg == '-ro'):
        normalize_by_ro = True
    

# Get basic grid info
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr = len(rr)

# radius normalized by r_o
rr_ro = rr/ro

# radius normalized by solar radius
rr_rsun = rr/6.96e10

ir_tbl = np.argmin(np.abs(rr_rsun - r_tbl))
if (normalize_by_ro):
    ir_tbl = np.argmin(np.abs(rr_ro - r_tbl))

np.save(datadir + dirname_stripped + '_ir_tbl.npy', ir_tbl)