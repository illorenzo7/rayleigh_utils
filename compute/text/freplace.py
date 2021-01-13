# Author: Loren Matilsky
# Created: 01/12/2021
# 
# usage: python $raco/freplace.py dirname .ext block1.txt block2.txt
# 
# searches [dirname] / files with .ext extension block1.txt
# and replace with block2.txt
#
# if -l is specified, will search line by line; if block1 is contained in
# line, will replace line with block2 
# (block1 and block2 should both be one line)

import sys, os
args = sys.argv[1:]
nargs = len(args)
dirname = args[0]
ext = args[1]
blockfile1 = args[2]
blockfile2 = args[3]

linebyline = False
for i in range(nargs):
    arg = args[i]
    if arg == '-l':
        linebyline = True

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
        if linebyline:
            block1 = block1.replace('\n','')
            st = open(dirname + '/' + fname).read()
            lines = open(dirname + '/' + fname).readlines()
            for line in lines:
                if block1 in line:
                    fnames_changed.append(dirname + '/' + fname)
                    stnew = st.replace(line, block2)
                    f = open(dirname + '/' + fname, "w")
                    f.write(stnew)
                    f.close()
                    
        else:
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
