# mark all subdirectories of desired directory 
# (argument 1, as a relative path to current working diretory)
# as backed up
# (put file 00_backedup in subdirectory)
#!/bin/bash

# Read in desired directory and get its subdirectories
maindir=$1
subdirs=`find $maindir -maxdepth 1 -mindepth 1 -type d -printf '%f\n'`

for subdir in $subdirs
do
    echo "touching $maindir/$subdir/00_backedup"
    touch $maindir/$subdir/00_backedup
done
# For good measure mark the main directory as all backedup
echo "touching $maindir/00_all_backedup"
touch $maindir/00_all_backedup
