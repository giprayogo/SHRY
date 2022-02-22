#!/bin/tcsh
#PBS -N prep
#PBS -j oe
#PBS -q SINGLE
#PBS -l select=1:ncpus=128:mpiprocs=1

source /etc/profile.d/modules.csh
module purge
module load oneapi-intel

setenv OMP_NUM_THREADS 1
setenv OPENBLAS_NUM_THREADS 1
setenv MKL_NUM_THREADS 1
setenv VECLIB_MAXIMUM_THREADS 1
setenv NUMEXPR_NUM_THREADS 1

cd ${PBS_O_WORKDIR}
#cat ${PBS_NODEFILE} > nodelist

#activate shry
setenv PATH "/home/nkousuke/application/anaconda3/bin:$PATH"
conda activate shry

set SCRIPT="."
set xsl_file="full_consistency_check.xls"

#prep_structure.py
python $SCRIPT/prep_structure.py $xsl_file > out_prep_structure

exit
