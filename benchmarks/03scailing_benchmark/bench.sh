#!/bin/tcsh
#PBS -N sg_sg_input
#PBS -j oe
#PBS -q SINGLE
#PBS -l select=1:ncpus=16:mpiprocs=1

source /etc/profile.d/modules.csh

module purge
module load PrgEnv-intel/2020
setenv OMP_NUM_THREADS 1

cd ${PBS_O_WORKDIR}
#cat ${PBS_NODEFILE} > nodelist

#activate shry
setenv PATH "/home/nkousuke/application/anaconda3/bin:$PATH"
conda activate shry

set sg=sg_input
set xsl_file="benchmark_SG_${sg}.xls"
python bench.py $xsl_file SG${sg} > out_bench_SG${sg}

exit
