#!/bin/bash
#PBS -N liq1
#PBS -q workq
#PBS -l select=4:ncpus=24:mpiprocs=24:interconnect=fdr:mem=60gb,walltime=72:00:00
#PBS -j oe

echo "START-------------------------"
qstat -xf $PBS_JOBID

module purge
module add gcc/7.1.0 openmpi/1.8.4 fftw/3.3.4-g481 

cd $PBS_O_WORKDIR

mpiexec -n 96 /home/xiaohoz/bin/lmp_mpi < input.prod


rm -f core.*



