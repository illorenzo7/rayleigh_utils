# Author: Loren Matilsky
# Created: 01/12/2021
# 
# usage: python $raco/freplace.py dirname .ext block1.txt block2.txt
# 
# searches [dirname] / files with .ext extension block1.txt
# and replace with block2.txt

import sys, os
dirname = sys.argv[1]
ext = sys.argv[2]
blockfile1 = sys.argv[3]
blockfile2 = sys.argv[4]

block1 = open(blockfile1).read()
block2 = open(blockfile2).read()

fnames_all = os.listdir(dirname)
fnames = []
for fname in fnames_all:
    if ext in fname:
        fnames.append(fname)

print ("=========================================")
print ("scanning %s files in %s" %(ext, dirname))
print ("=========================================")
print ("replacing")
print ("=========================================")
print (block1)
print ("with")
print ("=========================================")
print (block2)
print ("=========================================")

print ("cannot open:")
print ("-----------------------------------------")
fnames_changed = []
for fname in fnames:
    try:
        st = open(dirname + '/' + fname).read()
        if block1 in st:
            fnames_changed.append(dirname + '/' + fname)
            stnew = st.replace(block1, block2)
            f = open(dirname + '/' + fname, "w")
            f.write(stnew)
            f.close()
    except:
        print (dirname + '/' + fname)

print ("=========================================")
print ("found and replaced in:")
print ("-----------------------------------------")
for fname in fnames_changed:
    print(fname)
