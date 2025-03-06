# Author: Loren Matilsky
# Created: 01/12/2021
# 
# usage: python $ratx/freplace.py [dirname] [ext] [string1] [string2]
# can also replace string1 and/or string2 by b1 and/or b2
# which will use $ratx/block1.txt and $ratx/block2.txt instead of the strings
# 
# searches [dirname] for files ending in [ext] and replaces all occurrences of
# string1 (or b1) with string2 (or b2)

import sys, os
sys.path.append(os.environ['raco'])
from common import buff_line
args = sys.argv[1:]
nargs = len(args)
dirname = args[0]
ext = args[1]

string1 = args[2]
string2 = args[3]

if string1 == 'b1':
    string1 = open(os.environ['ratx'] + '/block1.txt').read()
if string2 == 'b2':
    string2 = open(os.environ['ratx'] + '/block2.txt').read()

# maybe do a "dry" run
dry = False
for arg in args:
    if arg == "--dry":
        dry = True
        print(buff_line)
        print("DOING A DRY RUN")

# loop through files, find and replace
fnames_all = []
for path, subdirs, files in os.walk(dirname):
    for name in files:
        fnames_all.append(os.path.join(path, name))
fnames = []
for fname in fnames_all:
    if ext == fname[-len(ext):]:
        fnames.append(fname)

print(buff_line)
print ("scanning %s files in %s" %(ext, dirname))
if dry:
    print ("would replace all occurrences of ")
else:
    print ("replacing all occurrences of ")
print (string1)
print ("with")
print (string2)

fnames_changed = []
for fname in fnames:
    try:
        st = open(fname).read()
        if string1 in st:
            fnames_changed.append(fname)
            if not dry:
                stnew = st.replace(string1, string2)
                f = open(fname, "w")
                f.write(stnew)
                f.close()
    except:
        print ("cannot read ", fname)

print(buff_line)
if dry:
    print ("would find and replace in the following files:")
else:
    print ("found and replaced in the following files:")

for fname in fnames_changed:
    print(fname)
