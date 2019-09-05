# Author: Loren Matilsky
# Created: 04/05/2019
# Takes sequential images from one directory (dir1) and copies every nth 
# one to another directory (dir2). 
# Usage: python $co/thin_images.py [dir1] [dir2] [n]
import sys, os
from shutil import copy

dir1, dir2, n = sys.argv[1:]
n = int(n)

all_files = os.listdir(dir1)
img_files = []
for fname in all_files:
    if '.png' in fname:
        img_files.append(fname)
img_files.sort()

for i in range(len(img_files)):
    fname = img_files[i]
    if i%n == 0:
        copy(dir1 + '/' + fname, dir2 + '/' + fname)
