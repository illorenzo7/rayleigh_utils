# Author: Loren Matilsky
# Created: well before 05/06/2019
# Takes the folder dirname (1st CLA) and renames .png files in sequential
# order, using the prefix 'img' and integers filled with 4 zeros, e.g., 
# img0000.png, img0001.png, etc.

import sys, os

dirname = sys.argv[1]

all_files = os.listdir(dirname)
img_files = []
for fname in all_files:
    if '.png' in fname:
        img_files.append(fname)

img_files.sort()
# First rename the files which may already have similar names to the new 
# ones (makes sure nothing gets moved twice, overwriting files)

count = 0
for fname in img_files:
    newname = 'prefix-' + fname
    os.rename(dirname + '/' + fname, dirname + '/' + newname)
    img_files[count] = newname
    count +=1

img_files.sort()
count = 0
for fname in img_files:
    newname = 'img%04i.png' %count
    os.rename(dirname + '/' + fname, dirname + '/' + newname)
    count += 1
