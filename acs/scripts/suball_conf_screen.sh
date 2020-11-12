#!/bin/bash

while (( "$#" )); do
	 sbatch -J $1 submit_conf_screen.sh $1
	 shift
 done
