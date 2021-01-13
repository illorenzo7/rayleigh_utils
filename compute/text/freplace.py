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

block1 = open(blockfile1).readlines()
nlines1 = len(block1)
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
for line in block1:
    print (line, end='')
print ("with")
print ("=========================================")
for line in block2:
    print (line, end='')
print ("=========================================")

print ("cannot open:")
print ("-----------------------------------------")

fnames_changed = []
for fname in fnames:
    try:
        st = open(dirname + '/' + fname).read()
        lines = open(dirname + '/' + fname).readlines()
        nlines = len(lines)
        for i in range(len(lines) - nlines1 + 1):
            lines_to_compare = lines[i:i+nlines1]
            the_truth = True
            for j in range(nlines1):
                the_truth *= block1[j][:-1] in lines_to_compare[j]
            if the_truth and lines_to_compare != block1:
                # replace the correponding block
                print ("lines_to_compare = ", lines_to_compare)
                print ("block1 = ", block1)
                fnames_changed.append(dirname + '/' + fname)
                block_we_compared = ''
                for j in range(nlines1):
                    block_we_compared += lines_to_compare[j]
                stnew = st.replace(block_we_compared, block2)
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
