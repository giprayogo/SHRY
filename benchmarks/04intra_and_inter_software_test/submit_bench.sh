#!/bin/tcsh
#PBS -N bench
#PBS -j oe
#PBS -q SINGLE
#PBS -l select=1:ncpus=128:mpiprocs=1

source /etc/profile.d/modules.csh
module purge
module load oneapi-intel

cd ${PBS_O_WORKDIR}
#cat ${PBS_NODEFILE} > nodelist

#activate shry
setenv PATH "/home/nkousuke/application/anaconda3/bin:$PATH"
conda activate shry

set xsl_file="full_consistency_check.xls"
python bench.py $xsl_file SG_all > out_bench

exit
