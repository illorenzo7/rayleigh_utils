echo "#PBS -S /bin/bash" > lscript
echo "#PBS -N $1" >> lscript

modeltype=$2
if [ $modeltype == 'bro_ele' ]
then
    ncpus=28
    execext='.avx2'
elif [ $modeltype == 'sky_ele' ] || [ $modeltype == 'cas_ait' ] 
then
    ncpus=40
    execext='.avx512'
elif [ $modeltype == 'rom_ait' ] 
then
    ncpus=128
    execext='.rome.2.28'
elif [ $modeltype == 'mil_ait' ] 
then
    ncpus=128
    execext='.milan.2.28'
else
    echo "unknown model type $modeltype"
fi

nprocs=$3

select=$(($nprocs/$ncpus + 1))
if [ $(($nprocs%$ncpus)) == 0 ]
then
    select=$(($select - 1))
fi

nprow=`python $rau/bash/nprow.py $nprocs`
npcol=`python $rau/bash/npcol.py $nprocs`

group=${4:-s7614} # charge it to s7614 by default

echo "#PBS -l select=$select:ncpus=$ncpus:model=$modeltype" >> lscript
echo "#PBS -q long" >> lscript
echo "#PBS -l walltime=120:00:00" >> lscript
echo "#PBS -j oe" >> lscript
echo "#PBS -W group_list=$group" >> lscript
echo "#PBS -m e" >> lscript
echo "#PBS -l site=needed=/home1+/nobackupp17" >> lscript

echo >> lscript

echo "module purge" >> lscript

if [ $modeltype == 'bro_ele' ] || [ $modeltype == 'sky_ele' ] || [ $modeltype == 'cas_ait' ] 
then
    echo "module load mpi-hpe" >> lscript
    echo "module load comp-intel" >> lscript
elif [ $modeltype == 'sky_ele' ] || [ $modeltype == 'sky_ele' ] 
then
    echo "module load gcc/13.2 mpi-hpe/mpt" >> lscript
    echo "fftwlib=/home4/nfeather/software/fftw-3.3.10-Rome/" >> lscript
    echo "blaslib=/home4/nfeather/software/OpenBLAS-0.3.27-Rome/" >> lscript
    echo "export LD_LIBRARY_PATH=$blaslib\lib:$fftwlib\lib:$LD_LIBRARY_PATH" >> lscript
    echo "export OMP_NUM_THREADS=1" >> lscript
fi

echo >> lscript

echo "export MPI_LAUNCH_TIMEOUT=40" >> lscript
echo "export PATH=\$PATH:/u/scicon/tools/bin" >> lscript
echo "export MPI_BUFS_PER_PROC=256" >> lscript
echo "export MPI_CHECK_ARGS=1" >> lscript

echo >> lscript

echo "cd \$PBS_O_WORKDIR" >> lscript

echo >> lscript

echo "/u/scicon/tools/bin/several_tries mpiexec -np $nprocs ./rayleigh$execext -nprow $nprow -npcol $npcol >> logfile2" >> lscript
