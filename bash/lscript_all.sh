# make 4 $fname files, one for each node type
declare -a arr=('bro_ele' 'sky_ele' 'cas_ait' 'rom_ait' 'mil_ait')

runname=$1
nprocs=$2
group=${3:-s7614} # charge it to s7614 by default

for modeltype in "${arr[@]}"
do 
    fname="lscript_$modeltype"
    bash /home1/lmatilsk/rayleigh/utils/bash/lscript.sh $fname $runname $modeltype $nprocs $group
done
