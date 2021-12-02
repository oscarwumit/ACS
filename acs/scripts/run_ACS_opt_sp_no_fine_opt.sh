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

step6="submit_molpro_array.sh"

step7="${ACS_PATH}acs/analysis/analyze_sp_result.py"

source activate acs_env

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


# step 6: run single point calculations
echo "Step 6: Submitting single point job..."
cd ../sp
START_TIME=$SECONDS
sbatch -W $step6
ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "Completed step 6 in $ELAPSED_TIME seconds" $'\n'


# step 7: analyze the single point calculations
echo "Step 7: Analyzing single point job..."
START_TIME=$SECONDS
python $create_generic_submit_script --step 9
sbatch -W sp_project_info.sh
ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo "Completed step 7 in $ELAPSED_TIME seconds" $'\n'

