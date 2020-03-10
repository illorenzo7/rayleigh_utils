# back up all subdirectories of desired directory 
# (argument 1, as a relative path to current working diretory)
# to lou (argument 2, as a full path)
# use shiftc and create .tar files
#!/bin/bash

# Read in desired directory and get its subdirectories
maindir=$1
fullpath=`readlink -f $maindir`
backupdir=$2
subdirs=`find $maindir -maxdepth 1 -mindepth 1 -type d -printf '%f\n'`

currentdir=`pwd`

# Make the corresponding directory and subdirectories on lou,
# also launch transfers in the loop
mkdir -p $backupdir
#fullpath=`readlink -f $maindir`
#remove="/noabackupp2/lmatilsk"
#backupdir=${fullpath//$remove/}
for subdir in $subdirs
do
    pwd
    mkdir $backupdir/$subdir
    cd $fullpath/$subdir
    shiftc --create-tar --index-tar --hosts=6 * $backupdir/$subdir/$subdir.tar
done
cd $currentdir
