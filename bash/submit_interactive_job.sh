thescript="/home1/lmatilsk/rayleigh/utils/bash/the_script_interactive"
echo "#PBS -S /bin/bash" > $thescript
echo "#PBS -N interactive_quote_unquote" >> $thescript

modeltype=$1
if [ $modeltype == 'bro_ele' ]
then
    ncpus=28
elif [ $modeltype == 'sky_ele' ] || [ $modeltype == 'cas_ait' ] 
then
    ncpus=40
elif [ $modeltype == 'rom_ait' ] 
then
    ncpus=128
elif [ $modeltype == 'mil_ait' ] 
then
    ncpus=128
else
    echo "unknown model type $modeltype"
fi

select=$2
thequeue=${3:-vlong}

echo "#PBS -l select=$select:ncpus=$ncpus:model=$modeltype" >> $thescript
echo "#PBS -q $thequeue" >> $thescript
if [[ "$thequeue" == "vlong" ]]; then
    echo "#PBS -l walltime=384:00:00" >> $thescript
elif [[ "$thequeue" == "long" ]]; then
    echo "#PBS -l walltime=120:00:00" >> $thescript
elif [[ "$thequeue" == "devel" ]]; then
    echo "#PBS -l walltime=2:00:00" >> $thescript
elif [[ "$thequeue" == "debug" ]]; then
    echo "#PBS -l walltime=2:00:00" >> $thescript
elif [[ "$thequeue" == "normal" ]]; then
    echo "#PBS -l walltime=4:00:00" >> $thescript
else
    echo "$thequeue is not a valid queue name. Aborting..."
    exit 1
fi

echo "#PBS -j oe" >> $thescript
echo "#PBS -W group_list=s7614" >> $thescript
echo "#PBS -m e" >> $thescript
echo "#PBS -l site=needed=/home1+/nobackupp17" >> $thescript

echo >> $thescript

echo "module purge" >> $thescript
echo "module load mpi-hpe" >> $thescript
echo "module load comp-intel" >> $thescript

echo >> $thescript
echo "cat \$PBS_NODEFILE > /nobackup/lmatilsk/nodefile_\$PBS_JOBID" >> $thescript
echo >> $thescript

echo "sleep 432000" >> $thescript

qsub $thescript
