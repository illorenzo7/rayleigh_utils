# extract all subdirectories of desired directory 
# (argument 1, as a relative path to current working diretory)
# from lou 
# to /nobackup subdir (argument 2, as a full path)
# use shiftc; write a temporary text file
# then extract .tar files
#!/bin/bash

# Read in desired directory to extract and get its subdirectories
extractdir=$1
subdirs=`find $extractdir -maxdepth 1 -mindepth 1 -type d -printf '%f\n'`

# get the directory to transfer to and create it
transferdir=$2
mkdir -p $nb/$transferdir

# make a file of transfer instructions (a list of sources and dest)
thefile=~/transfer_instructions_tmp
echo -n '' > $thefile # erase the current file
# -n means don't add newline character, which is added by echo by default

# loop over the subdirs and add shift sources and destinations
# one line at a time
for subdir in $subdirs
do
    echo -n $extractdir/$subdir/\*.tar >> $thefile # extraction dir
    echo -n ' ' >> $thefile # space
    echo $nb/$transferdir/$subdir/ >> $thefile # transfer directory
    mkdir $nb/$transferdir/$subdir # make the transfer directory
done

# now the actual transfer
echo $bufferstring
echo "about to run shiftc --extract-tar --hosts=6"
echo "on $thefile, which is"
cat $thefile
echo $bufferstring
shiftc --extract-tar --hosts=6 < $thefile
echo $bufferstring
