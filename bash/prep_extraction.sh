# Author: Loren Matilsky
# Created: 12/04/2020
# find all subdirectories of maindir (argument 1) and create a text file
# subdirs.txt in maindir with one line for each subdirectory
# this will allow remote extraction of all subdirectories using multiple 
# shiftc transfers
#!/bin/bash

# Read in desired directory and output its subdirectories to a .txt file
maindir=$1
fullpath=`readlink -f $maindir`
find $maindir -maxdepth 1 -mindepth 1 -type d -printf '%f\n' > $fullpath/subdirs.txt
