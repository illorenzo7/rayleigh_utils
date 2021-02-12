# Loren Matilsky
# 07/06/2018

# Takes "dirname" as argument and for all files in "dirname", replaces spaces
# in the filename with underscores

import os, sys

dry = False
dirname = sys.argv[1]
other_args = sys.argv[2:]
nother_args = len(other_args)
for i in range(nother_args):
    if other_args[i] == '-dry':
        dry = True

all_files = os.listdir(dirname)


print ('In %s:' %dirname)
for fname in all_files:
    new_fname = ''
    nchar = len(fname)
    for i in range(nchar):
        current_char = fname[i]
        
        if not (current_char == ' ' and fname[i+1:] == '.pdf'):
            if current_char == ' ':
                new_fname += '_'
            else:
                new_fname += current_char

    
    if (new_fname != fname):
        if not dry:
            os.rename(dirname + '/' + fname, dirname + '/' + new_fname)
            print ('Moving %s to\n %s...' %(fname, new_fname))
            print ('-----------------------------')
        else:
            print ('Will move %s to\n %s...' %(fname, new_fname))
            print ('-----------------------------')            