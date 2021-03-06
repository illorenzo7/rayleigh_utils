# count files and subdirectories of all subdirectories of 
# desired directory using "tree"
# (argument 1, as a relative path to current working diretory)

# Read in desired directory and get its subdirectories
maindir=$1
fullpath=`readlink -f $maindir`
subdirs=`find $maindir -maxdepth 1 -mindepth 1 -type d -printf '%f\n'`

# Return to our current location at the end of the script
currentdir=`pwd`

echo "===================================================================="
echo "===================================================================="
echo "In $fullpath/:"
tree -a . | tail -n 1
echo "===================================================================="
echo "===================================================================="

# First get the number of files in the whole directory

# Loop over subdirectories and count their contents with "tree"
for subdir in $subdirs
do
    cd $fullpath/$subdir
    echo "In subdir $subdir/:"
    tree -a . | tail -n 1
    echo "==============================="
done
cd $currentdir
