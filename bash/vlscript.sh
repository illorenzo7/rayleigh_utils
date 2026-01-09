fname=$1
echo "#PBS -S /bin/bash" > $fname
echo "#PBS -N $2" >> $fname

modeltype=$3
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
    execext='.rome2.28'
elif [ $modeltype == 'mil_ait' ] 
then
    ncpus=128
    execext='.milan2.30'
else
    echo "unknown model type $modeltype"
fi

nprocs=$4

select=$(($nprocs/$ncpus + 1))
if [ $(($nprocs%$ncpus)) == 0 ]
then
    select=$(($select - 1))
fi

nprow=`python $rau/bash/nprow.py $nprocs`
npcol=`python $rau/bash/npcol.py $nprocs`

group=${5:-s7614} # charge it to s7614 by default

echo "#PBS -l select=$select:ncpus=$ncpus:model=$modeltype" >> $fname
echo "#PBS -q vlong" >> $fname
echo "#PBS -l walltime=384:00:00" >> $fname
echo "#PBS -j oe" >> $fname
echo "#PBS -W group_list=$group" >> $fname
echo "#PBS -m e" >> $fname
echo "#PBS -l site=needed=/home1+/nobackupp17" >> $fname

echo >> $fname

echo "module purge" >> $fname

if [ $modeltype == 'bro_ele' ] || [ $modeltype == 'sky_ele' ] || [ $modeltype == 'cas_ait' ] 
then
    echo "module load mpi-hpe" >> $fname
    echo "module load comp-intel" >> $fname
elif [ $modeltype == 'rom_ait' ] 
then
    echo "module load gcc/13.2 mpi-hpe/mpt.2.28_25Apr23_rhel87" >> $fname
elif [ $modeltype == 'mil_ait' ] 
then
    echo "module load gcc/13.2 mpi-hpe/mpt" >> $fname
fi

if [ $modeltype == 'rom_ait' ] || [ $modeltype == 'mil_ait' ] 
then
    echo "fftwlib=/home4/nfeather/software/fftw-3.3.10-Rome/" >> $fname
    echo "blaslib=/home4/nfeather/software/OpenBLAS-0.3.27-Rome/" >> $fname
    echo 'export LD_LIBRARY_PATH=$blaslib\lib:$fftwlib\lib:$LD_LIBRARY_PATH' >> $fname
    echo "export OMP_NUM_THREADS=1" >> $fname
fi

echo >> $fname

echo "export MPI_LAUNCH_TIMEOUT=40" >> $fname
echo "export PATH=\$PATH:/u/scicon/tools/bin" >> $fname
echo "export MPI_BUFS_PER_PROC=256" >> $fname
echo "export MPI_CHECK_ARGS=1" >> $fname

echo >> $fname

echo "cd \$PBS_O_WORKDIR" >> $fname

echo >> $fname

echo "/u/scicon/tools/bin/several_tries mpiexec -np $nprocs ./rayleigh$execext -nprow $nprow -npcol $npcol >> logfile2" >> $fname
