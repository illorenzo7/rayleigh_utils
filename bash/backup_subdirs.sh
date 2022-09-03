# back up all subdirectories of desired directory 
# (argument 1, as a relative path to current working diretory)
# to lou (argument 2, as a full path)
# use shiftc and create .tar files
#!/bin/bash

# Read in desired directory and get its subdirectories
maindir=$1
#fullpath=`readlink -f $maindir`
#cd $fullpath
subdirs=`find $maindir -maxdepth 1 -mindepth 1 -type d -printf '%f\n'`

# move in to the desired directory
cd $maindir
echo $bufferstring
echo -n "I am here: "
echo `pwd`
# remember current directory to cd back into
currentdir=`pwd`

# get the backup directory and create it
backupdir=$2
mkdir -p $backupdir

# Make the corresponding directory and subdirectories on lou,
## also launch transfers in the loop
#fullpath=`readlink -f $maindir`
#remove="/noabackupp2/lmatilsk"
#backupdir=${fullpath//$remove/}

# make a file of transfer instructions (a list of sources and dest)
thefile=~/transfer_instructions_tmp
echo -n '' > $thefile
for subdir in $subdirs
do
    #pwd
    #mkdir $backupdir/$subdir
    #cd $fullpath/$subdir
    #ls | tr "\n" " " >> ~/transfer_instructions_tmp
    #shiftc --create-tar --index-tar --hosts=6 * $backupdir/$subdir/$subdir.tar
    echo -n $subdir >> $thefile
    echo -n ' ' >> $thefile
    echo $backupdir/$subdir.tar >> $thefile
done

# now the actual transfer
echo $bufferstring
echo "running shiftc --create-tar --index-tar --hosts=6"
echo "on $thefile, which is"
cat $thefile
echo $bufferstring
shiftc --create-tar --index-tar --hosts=6 $thefile
cd $currentdir
echo $bufferstring
echo -n "I returned here: "
echo `pwd`
echo $bufferstring
