#!/bin/bash

# my "mirror" directory
maindir=`pwd`

# get all subdirs
cd $maindir
subdirs=`ls -d Case*`

source /home1/lmatilsk/environment/dotbash_profile
source /home1/lmatilsk/environment/dotbashrc
source /home1/lmatilsk/environment/bashrc_lfe

# copy over the subdirs and softlink their contents
for subdir in $subdirs
do
    #echo "cd $maindir"
    cd $maindir/$subdir
    pwd
    mpiexec -n 72 python -u /home1/lmatilsk/rayleigh/utils/compute/timetrace/G_Avgs_quad.py . --all --thinby 50
    #mpiexec -n 72 python -u /home1/lmatilsk/rayleigh/utils/zz_legacy/raslice_nondim_specificsuite.py . --varnames vrprime omzprime bp --samplevals all
    #python /home1/lmatilsk/rayleigh/utils/plot/timetrace/etrace_quad.py . --noshow --log
    #python /home1/lmatilsk/rayleigh/utils/zz_legacy/mercirc_nondim_specificsuite.py . --rvals 5e10 --noshow
    #python /home1/lmatilsk/rayleigh/utils/plot/azav/diffrot.py . --noshow --rvals 5e10

done

# return home
cd $maindir
