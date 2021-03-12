#!/bin/bash
#SBATCH -J psi4
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --time=6:00:00

ARCPATH=/home/gridsan/kspieker/RMG/ARC
RMGPATH=/home/gridsan/kspieker/RMG/RMG-Py
ACSPATH=/home/gridsan/kspieker/RMG/ACS
export PYTHONPATH=$PYTHONPATH:$RMGPATH
export PYTHONPATH=$PYTHONPATH:$ARCPATH
export PYTHONPATH=$PYTHONPATH:$ACSPATH
export PSI_SCRATCH=/home/gridsan/kspieker/scratch/psi4


source activate arc_env

python /home/gridsan/kspieker/RMG/ACS/acs/screening/run_screening.py  $1
