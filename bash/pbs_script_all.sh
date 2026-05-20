# make 4 $fname files, one for each node type
declare -a arr=('bro_ele' 'sky_ele' 'cas_ait' 'rom_ait' 'mil_ait')

runname=$1
nprocs=$2
qtype=${3:-vlong} # create the vlong scripts by default
group=${4:-s3058} # charge it here by default
if [ $qtype == 'vlong' ]
then 
    basename="vlscript"
elif [ $qtype == 'long' ]
then
    basename="lscript"
elif [ $qtype == 'normal' ]
then
    basename="nscript"
elif [ $qtype == 'devel' ] 
then
    basename="dscript"
elif [ $qtype == 'debug' ] 
then
    basename="dbscript"
fi

for modeltype in "${arr[@]}"
do 
    fname=${basename}_${modeltype}
    bash /home1/lmatilsk/rayleigh/utils/bash/pbs_script.sh $fname $runname $modeltype $nprocs $qtype $group
done
