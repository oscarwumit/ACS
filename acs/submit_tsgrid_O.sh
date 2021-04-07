#!/bin/bash
#SBATCH -J psi4
#SBATCH -N 1
#SBATCH --exclusive

ARCPATH=~/Software/ARC
RMGPATH=~/Software/RMG-Py
export PYTHONPATH=$PYTHONPATH:$RMGPATH

source activate ACS_env

python /home/gridsan/oscarwu/GRPAPI/imipramine_site_5_2d_grid_sp/TS_confgrid.py  $1

source deactivate


