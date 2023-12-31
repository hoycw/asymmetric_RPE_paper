#!/bin/sh
# This script analyzes subject data in ~/MATLAB/DATA/Avgusta for each task
# The SGE_TASK variable specifies the data sets, as in $HOME/test.qsub

cd /home/knight/hoycw/PRJ_Error/scripts/

# define where the datasets are located
#PRCSDATADIR="/home/knight/hoycw/PRJ_Error/data/"

# don't change this variable
# used by the submit script to define which data sets to analyze
SBJ="${SGE_TASK}"

# define function
FUNCTION='SBJ02_preproc'

# set up matlab function call
func_call="${FUNCTION}('${SBJ}', '${proc_id}')"

# define commands to execute via SGE
echo ${SBJ}
echo ${func_call}
echo $$
echo ${func_call} > NotBackedUp/tmpSGE/${FUNCTION}_${SBJ}.m
time matlab -nodesktop -nosplash -nodisplay < NotBackedUp/tmpSGE/${FUNCTION}_${SBJ}.m
rm NotBackedUp/tmpSGE/${FUNCTION}_${SBJ}.m
