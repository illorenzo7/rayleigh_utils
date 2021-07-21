echo "#PBS -S /bin/bash" > lscript
echo "#PBS -N $1" >> lscript

modeltype=$2
if [ $modeltype == 'bro' ]
then
    ncpus=28
elif [ $modeltype == 'has' ]
then
    ncpus=24
elif [ $modeltype == 'ivy' ] 
then
    ncpus=20
elif [ $modeltype == 'san' ] 
then
    ncpus=16
else
    echo "unknown model type $modeltype"
fi

nprocs=$3

select=$(($nprocs/$ncpus + 1))
if [ $nprocs%$ncpus == 0 ]
then
    select=$(($select + 1))
fi

nprow=`python $rau/zz_runscripts/nprow.py $nprocs`
npcol=`python $rau/zz_runscripts/npcol.py $nprocs`

echo "#PBS -l select=$select:ncpus=$ncpus:model=$modeltype" >> lscript
echo "#PBS -q long" >> lscript
echo "#PBS -l walltime=120:00:00" >> lscript
echo "#PBS -j oe" >> lscript
echo "#PBS -W group_list=s2051" >> lscript
echo "#PBS -m e" >> lscript
echo "#PBS -l site=needed=/home1+/nobackupp18" >> lscript

echo >> lscript

echo "module purge" >> lscript
echo "module load mpi-hpe" >> lscript
echo "module load comp-intel" >> lscript

echo >> lscript

echo "export MPI_LAUNCH_TIMEOUT=40" >> lscript
echo "export PATH=\$PATH:/u/scicon/tools/bin" >> lscript
echo "export MPI_BUFS_PER_PROC=256" >> lscript
echo "export MPI_CHECK_ARGS=1" >> lscript

echo >> lscript

echo "cd \$PBS_O_WORKDIR" >> lscript

echo >> lscript

echo "/u/scicon/tools/bin/several_tries mpiexec -np 1024 ./rayleigh.opt -nprow $nprow -npcol $npcol >> logfile2" >> lscript
