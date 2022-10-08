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
string1 = args[2]
string2 = args[3]

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
print ("replacing all occurrences of ")
print ("=========================================")
print (string1)
print ("with")
print ("=========================================")
print (string2)
print ("=========================================")

print ("cannot open:")
print ("-----------------------------------------")

fnames_changed = []
for fname in fnames:
    try:
        st = open(fname).read()
        if string1 in st:
            fnames_changed.append(fname)
            stnew = st.replace(string1, string2)
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
