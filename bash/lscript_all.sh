# make 4 $fname files, one for each node type
declare -a arr=('bro' 'has' 'ivy' 'san' 'sky_ele' 'cas_ait' 'mil_ait' 'rom_ait')

for modeltype in "${arr[@]}"
do 
    fname="lscript_$modeltype"
    echo "#PBS -S /bin/bash" > $fname
    echo "#PBS -N $1" >> $fname

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
    elif [ $modeltype == 'mil_ait' ] || [ $modeltype == 'rom_ait' ] 
    then
        ncpus=128
        execext='.rome2.28'
    fi

    nprocs=$2

    select=$(($nprocs/$ncpus + 1))
    if [ $(($nprocs%$ncpus)) == 0 ]
    then
        select=$(($select - 1))
    fi

    nprow=`python $rau/bash/nprow.py $nprocs`
    npcol=`python $rau/bash/npcol.py $nprocs`

    echo "#PBS -l select=$select:ncpus=$ncpus:model=$modeltype" >> $fname
    echo "#PBS -q long" >> $fname
    echo "#PBS -l walltime=120:00:00" >> $fname
    echo "#PBS -j oe" >> $fname
    echo "#PBS -W group_list=s2328" >> $fname
    echo "#PBS -m e" >> $fname
    echo "#PBS -l site=needed=/home1+/nobackupp17" >> $fname

    echo >> $fname

    echo "module purge" >> $fname
    if [ $modeltype == 'mil_ait' ] || [ $modeltype == 'rom_ait' ]
    then
        echo "module load gcc/13.2 mpi-hpe/mpt" >> $fname
        echo "fftwlib=/home4/nfeather/software/fftw-3.3.10-Rome/" >> $fname
        echo "blaslib=/home4/nfeather/software/OpenBLAS-0.3.27-Rome/" >> $fname
        echo "export LD_LIBRARY_PATH=\$blaslib\\lib:\$fftwlib\\lib:\$LD_LIBRARY_PATH " >> $fname
    else
        echo "module load mpi-hpe" >> $fname
        echo "module load comp-intel" >> $fname
    fi

    echo >> $fname

    echo "export OMP_NUM_THREADS=1" >> $fname
    echo "export MPI_LAUNCH_TIMEOUT=40" >> $fname
    echo "export PATH=\$PATH:/u/scicon/tools/bin" >> $fname
    echo "export MPI_BUFS_PER_PROC=256" >> $fname
    echo "export MPI_CHECK_ARGS=1" >> $fname

    echo >> $fname

    echo "cd \$PBS_O_WORKDIR" >> $fname

    echo >> $fname

    echo "/u/scicon/tools/bin/several_tries mpiexec -np $nprocs ./rayleigh$execext -nprow $nprow -npcol $npcol >> logfile2" >> $fname
done
