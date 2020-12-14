#!/bin/bash

ACS_PATH="/home/gridsan/kspieker/ENI/ACS"

while (( "$#" )); do
	 sbatch -W -J $1 "${ACS_PATH}/acs/scripts/submit_conf_screen.sh" $1
	 shift
 done
wait
