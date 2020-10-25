# remove all but last 5 Shell_Slices and Shell_Spectra
# from [dirname] (argument 1, as a relative path to current working diretory)
#!/bin/bash

# Read in desired directory and get current directory
maindir=$1
odir=`pwd`

# Remove the Shell_Slices
nfiles=`ls $maindir/Shell_Slices | wc -l`
allbut5=`ls $maindir/Shell_Slices | head -n $((nfiles -5))`
cd $maindir/Shell_Slices
rm $allbut5
cd $odir

# Remove the Shell_Slices
nfiles=`ls $maindir/Shell_Spectra | wc -l`
allbut5=`ls $maindir/Shell_Spectra | head -n $((nfiles -5))`
cd $maindir/Shell_Spectra
rm $allbut5
cd $odir
