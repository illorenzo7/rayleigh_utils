lscript="/home1/lmatilsk/rayleigh/utils/bash/lscript_interactive"
echo "#PBS -S /bin/bash" > $lscript
echo "#PBS -N interactive_quote_unquote" >> $lscript

modeltype=$1
if [ $modeltype == 'bro' ]
then
    ncpus=28
    execext='.avx2'
elif [ $modeltype == 'has' ]
then
    ncpus=24
    execext='.avx2'
elif [ $modeltype == 'ivy' ] 
then
    ncpus=20
    execext='.avx'
elif [ $modeltype == 'san' ] 
then
    ncpus=16
    execext='.avx'
elif [ $modeltype == 'sky_ele' ] || [ $modeltype == 'cas_ait' ] 
then
    ncpus=40
    execext='.avx512'
else
    echo "unknown model type $modeltype"
fi

select=$2

echo "#PBS -l select=$select:ncpus=$ncpus:model=$modeltype" >> $lscript
echo "#PBS -q long" >> $lscript
echo "#PBS -l walltime=120:00:00" >> $lscript
echo "#PBS -j oe" >> $lscript
echo "#PBS -W group_list=s2328" >> $lscript
echo "#PBS -m e" >> $lscript
echo "#PBS -l site=needed=/home1+/nobackupp17" >> $lscript

echo >> $lscript

echo "module purge" >> $lscript
echo "module load mpi-hpe" >> $lscript
echo "module load comp-intel" >> $lscript

echo >> $lscript
echo "cat \$PBS_NODEFILE > /nobackup/lmatilsk/nodefile_\$PBS_JOBID" >> $lscript
echo >> $lscript

echo "sleep 432000" >> $lscript

qsub $lscript
