# make 4 $fname files, one for each node type
declare -a arr=('bro_ele' 'sky_ele' 'cas_ait' 'rom_ait' 'mil_ait')

nprocs=$2

for modeltype in "${arr[@]}"
do 
    fname="lscript_$modeltype"
    bash /home1/lmatilsk/rayleigh/utils/bash/lscript.sh $fname $modeltype $nprocs
done
