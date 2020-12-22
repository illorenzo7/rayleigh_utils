# Author: Loren Matilsky
# Created: 12/04/2020
# find all subdirectories of maindir (argument 1) on an lfe and extract all
# the data from the tar files to the backup dir (argument 2)
#!/bin/bash

# Read in desired (remote) directory and back up directory 
maindir=$1
extractdir=$2
# get the remote directory's subdirectories
#scp lfe:$maindir/ ls" > subdirs
# don't know why I can't scp the subdirs file from lou....

while read subdir; do
    mkdir $extractdir/$subdir
    echo "$extracting lfe:$maindir/$subdir/*.tar"
    sup shiftc --extract-tar lfe:$maindir/$subdir/*.tar $extractdir/$subdir/
done <$extractdir/subdirs.txt
