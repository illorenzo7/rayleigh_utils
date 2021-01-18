# Author: Loren Matilsky
# Created: 01/12/2021
# 
# usage: python $raco/freplace.py dirname .ext block1.txt block2.txt
# 
# searches [dirname] / files with .ext extension block1.txt
# and replace with block2.txt

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

fnames_all = []
for path, subdirs, files in os.walk(dirname):
    for name in files:
        fnames_all.append(os.path.join(path, name))
fnames = []
for fname in fnames_all:
    if ext == fname[-len(ext):]:
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
        st = open(fname).read()
        lines = open(fname).readlines()
        nlines = len(lines)
        for i in range(nlines - nlines1 + 1):
            lines_to_compare = lines[i:i+nlines1]
            the_truth = True
            for j in range(nlines1):
                the_truth *= block1[j][:-1] ==\
                    lines_to_compare[j][:len(block1[j]) - 1]

            if the_truth:
                # replace the correponding block if block 2 is different
                block_we_compared = ''
                for j in range(nlines1):
                    block_we_compared += lines_to_compare[j]
                if block2 != block_we_compared:
                    fnames_changed.append(fname)
                    stnew = st.replace(block_we_compared, block2)
                    f = open(fname, "w")
                    f.write(stnew)
                    f.close()
    except:
        print (fname)

print ("=========================================")
print ("found and replaced in:")
print ("-----------------------------------------")
for fname in fnames_changed:
    print(fname)
