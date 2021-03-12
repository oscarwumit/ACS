#!/bin/bash
#SBATCH -J run_ACS					# Job name
#SBATCH -o run_ACS.%N.%j.out      	# STDOUT
#SBATCH -e run_ACS.%N.%j.err      	# STDERR
#SBATCH -t 1-01:00:00 				# time (DD-HH:MM:SS)
#SBATCH -n 4          				# number of cores, default is one task per node
#SBATCH -N 1            			# number of nodes
#SBATCH --mem-per-cpu=4096        	# memory requested per cpu, default 2048 (optional)

# Activate Python environment and run program.

START_TIME=$SECONDS

# Run job
echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"


bash /home/gridsan/kspieker/RMG/ACS/acs/scripts/run_ACS.sh


ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "Elapsed time (s):" 
echo $ELAPSED_TIME

