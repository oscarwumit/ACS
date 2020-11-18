# step 0: define paths
##############################################
# only need to change this path
ACS_PATH="/home/gridsan/kspieker/ENI/ACS/"

##############################################

step1="${ACS_PATH}acs/screening/gen_conf.py"

step2="${ACS_PATH}acs/scripts/suball_conf_screen.sh"

step3="${ACS_PATH}acs/screening/analyze_screen.py"

step4="submit_qchem53_array.sh"

step5="${ACS_PATH}acs/analysis/analyze_opt_freq_result.py"

step6="submit_qchem53_array.sh"

step7="${ACS_PATH}acs/analysis/analyze_opt_freq_result.py"

step8="submit_orca_array.sh"

step9="${ACS_PATH}acs/analysis/analyze_sp_result.py"

step10="${ACS_PATH}acs/analysis/calc_mstst_rate.py"


# step 1: generate conformers
echo "Step 1: Generating conformers..."
source activate arc_env
python $step1 input.yml
echo "Completed step 1"


# step 2: run sp calculations with psi4
echo "Step 2: Running single point calculations with Psi4..."
cd initial_sp_screening
bash $step2 */*.yml
echo "Completed step 2"


# step 3: analyze sp screening results and find set of lowest energy conformers
# cd into one of the directories, such as 0, which will always be present
echo "Step 3: Analyzing Psi4 single point screening"
cd 0
python $step3 initial_conf_screening_result.yml
echo "Completed step 3"


# step 4: run geometry optimizations for the lower energy conformers
echo "Step 4: Running QM jobs..."
cd ../../opt
sbatch $step4
echo "Completed step 4"


# step 5: analyze geometry optimizations and frequency calculations
echo "Step 5: Analyzing results from QM job..."
python $step5 opt_project_info.yml
echo "Completed step 5"


# step 6: fine tune the optimization at higher level of theory
echo "Step 6: Running QM jobs for fine optimization..."
cd ../opt_fine
sbatch $step6
echo "Completed step 6"


# step 7: analzye the fine optimization
echo "Step 7: Analyzing results from fine QM job..."
python $step7 fine_opt_project_info.yml
echo "Completed step 7"


# step 8: run single point calculations
echo "Step 8: Submitting single point job..."
cd ../sp
bash $step8
echo "Completed step 8"


# step 9: analyze the single point calculations
echo "Step 9: Analyzing single point job..."
python $step9 sp_project_info.yml
echo "Completed step 9"


# step 10: run MSTST calculation with final set of optimized conformers
# this step needs the paths to r1, r2, p1, p2. not yet automated
# python step10

