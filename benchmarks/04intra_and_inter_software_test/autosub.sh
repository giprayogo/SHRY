#!/bin/tcsh
#PBS -N sg_sg_input
#PBS -j oe
#PBS -q DEFAULT
#PBS -l select=1:ncpus=64:mpiprocs=1

source /etc/profile.d/modules.csh
module purge
module load oneapi-intel/2021.1.1
module load oneapi/mpi/2021.1.1

setenv OMP_NUM_THREADS 1

cd ${PBS_O_WORKDIR}
#cat ${PBS_NODEFILE} > nodelist

#activate shry
setenv PATH "/home/nkousuke/application/anaconda3/bin:$PATH"
conda activate shry

#autosub
set sg=sg_input
python ./autosub.py SG${sg}> out_autosub_SG${sg}

