#!/bin/sh
# This script analyzes subject data in ~/MATLAB/DATA/Avgusta for each task
# The SGE_TASK variable specifies the data sets, as in $HOME/test.qsub

cd /home/knight/hoycw/PRJ_Error/scripts/

# define where the datasets are located
#PRCSDATADIR="/home/knight/hoycw/PRJ_Error/data/"

# don't change this variable
# used by the submit script to define which data sets to analyze
DATASET="${SGE_TASK}"

# define function
FUNCTION='SBJ00_cleaning_prep'

# set up matlab function call
func_call="${FUNCTION}('${DATASET}', '${plot_psd}')"

# define commands to execute via SGE
echo ${DATASET}
echo ${func_call}
echo $$
echo ${func_call} > NotBackedUp/tmpSGE/SBJ00_cleaning_prep_${DATASET}.m
time matlab -nodesktop -nosplash -nodisplay < NotBackedUp/tmpSGE/SBJ00_cleaning_prep_${DATASET}.m
rm NotBackedUp/tmpSGE/SBJ00_cleaning_prep_${DATASET}.m
