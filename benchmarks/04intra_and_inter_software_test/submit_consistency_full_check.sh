#!/bin/tcsh
#PBS -N full_SHRY
#PBS -j oe
#PBS -q XLARGE
#PBS -l select=24:ncpus=128:mpiprocs=64

source /etc/profile.d/modules.csh
module purge
module load oneapi-intel

setenv OMP_NUM_THREADS 1
setenv OPENBLAS_NUM_THREADS 1
setenv MKL_NUM_THREADS 1
setenv VECLIB_MAXIMUM_THREADS 1
setenv NUMEXPR_NUM_THREADS 1

cd ${PBS_O_WORKDIR}
cat ${PBS_NODEFILE} > nodelist

#activate shry
setenv PATH "/home/nkousuke/application/anaconda3/bin:$PATH"
conda activate shry

set SCRIPT="."
set xsl_file="full_consistency_check.xls"

#run check_full
mpirun -np 1536 --hostfile ${PBS_NODEFILE} python ${SCRIPT}/check_full_mpi.py $xsl_file > out_check_full_mpi

exit
