#!/bin/bash
#SBATCH -J psi4
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --time=6:00:00

ARCPATH=/home/gridsan/yunsie/git_repo/ARC
RMGPATH=/home/gridsan/yunsie/git_repo/RMG-Py
ACSPATH=/home/gridsan/yunsie/git_repo/ACS
export PYTHONPATH=$PYTHONPATH:$RMGPATH
export PYTHONPATH=$PYTHONPATH:$ARCPATH
export PYTHONPATH=$PYTHONPATH:$ACSPATH
export PSI_SCRATCH=/home/gridsan/yunsie/scratch/psi4


source activate ACS_env

python /home/gridsan/yunsie/git_repo/ACS/acs/screening/run_screening.py  $1
