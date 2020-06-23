#!/bin/bash

while (( "$#" )); do
	 sbatch -J $1 /home/gridsan/oscarwu/GRPAPI/imipramine_site_5_2d_grid_sp/submit_tsgrid_O.sh $1
	 shift
 done
