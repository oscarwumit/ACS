# step 0: define paths
##############################################
# only need to change this path
ACS_PATH="/home/gridsan/kspieker/RMG/ACS/"

##############################################

create_generic_submit_script="${ACS_PATH}acs/scripts/create_generic_submit_script.py"

step1="${ACS_PATH}acs/screening/gen_conf.py"

step2="submit_psi4_array.sh"

step3="${ACS_PATH}acs/screening/analyze_screen.py"

step4="submit_qchem53_array.sh"

step5="${ACS_PATH}acs/analysis/analyze_opt_freq_result.py"

step6="submit_qchem53_array.sh"

step7="${ACS_PATH}acs/analysis/analyze_opt_freq_result.py"

step8="submit_orca_array.sh"

step9="${ACS_PATH}acs/analysis/analyze_sp_result.py"

step10="${ACS_PATH}acs/analysis/calc_mstst_rate.py"


source activate arc_env

# step 1: generate conformers
echo "Step 1: Generating conformers..."
START_TIME=$SECONDS
python $create_generic_submit_script --step 1
sbatch -W gen_conf.sh
ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "Completed step 1 in $ELAPSED_TIME seconds" $'\n'


# step 2: run sp calculations with psi4
echo "Step 2: Running single point calculations with Psi4..."
cd initial_sp_screening
START_TIME=$SECONDS
python $create_generic_submit_script --step 2
sbatch -W $step2
ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "Completed step 2 in $ELAPSED_TIME seconds" $'\n'


# step 3: analyze sp screening results and find set of lowest energy conformers
# cd into one of the directories, such as 0, which will always be present
echo "Step 3: Analyzing Psi4 single point screening"
cd 0
START_TIME=$SECONDS
python $create_generic_submit_script --step 3
sbatch -W analyze_screen.sh
ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "Completed step 3 in $ELAPSED_TIME seconds" $'\n'


# step 4: run geometry optimizations for the lower energy conformers
echo "Step 4: Running QM jobs..."
cd ../../opt
START_TIME=$SECONDS
sbatch -W $step4
ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "Completed step 4 in $ELAPSED_TIME seconds" $'\n'


# step 5: analyze geometry optimizations and frequency calculations
echo "Step 5: Analyzing results from QM job..."
START_TIME=$SECONDS
python $create_generic_submit_script --step 5
sbatch -W opt_project_info.sh
ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "Completed step 5 in $ELAPSED_TIME seconds" $'\n'


# step 6: fine tune the optimization at higher level of theory
echo "Step 6: Running QM jobs for fine optimization..."
cd ../opt_fine
START_TIME=$SECONDS
sbatch -W $step6
ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "Completed step 6 in $ELAPSED_TIME seconds" $'\n'


# step 7: analzye the fine optimization
echo "Step 7: Analyzing results from fine QM job..."
START_TIME=$SECONDS
python $create_generic_submit_script --step 7
sbatch -W fine_opt_project_info.sh
ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "Completed step 7 in $ELAPSED_TIME seconds" $'\n'


# step 8: run single point calculations
echo "Step 8: Submitting single point job..."
cd ../sp
START_TIME=$SECONDS
sbatch -W $step8
ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "Completed step 8 in $ELAPSED_TIME seconds" $'\n'


# step 9: analyze the single point calculations
echo "Step 9: Analyzing single point job..."
START_TIME=$SECONDS
python $create_generic_submit_script --step 9
sbatch -W sp_project_info.sh
ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "Completed step 9 in $ELAPSED_TIME seconds" $'\n'


# step 10: run MSTST calculation with final set of optimized conformers
# this step needs the paths to r1, r2, p1, p2. not yet automated
# python step10

