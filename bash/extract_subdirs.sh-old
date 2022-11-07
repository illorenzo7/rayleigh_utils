# extract all subdirectories of desired directory on lou
# (argument 1, as full or relative path)
# to /nobackup (argument 2, as full or relative path)
# use shiftc and extract .tar files
#!/bin/bash

# Read in desired directory and get its subdirectories
extractdir=$1
extractdir=`readlink -f $extractdir`
maindir=$2 # directory put files (on /nobackup)
maindir=`readlink -f $maindir`
subdirs=`find $extractdir -maxdepth 1 -mindepth 1 -type d -printf '%f\n'`

# Make the extract directory's subdirectories on /nobackup
# also launch transfers in the loop

for subdir in $subdirs
do
    pwd
    mkdir $maindir/$subdir
    echo "shiftc --extract-tar --hosts=6 $extractdir/$subdir/*.tar $maindir/$subdir"
    shiftc --extract-tar --hosts=6 $extractdir/$subdir/*.tar $maindir/$subdir
    echo "======================================="
done
