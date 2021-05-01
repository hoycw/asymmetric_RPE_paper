#!/bin/sh

# Move to scripts directory 
root_dir="/Volumes/hoycw_clust/PRJ_Error/"
cd $root_dir/scripts/

# Load the correct SBJ list
echo "Please enter the desired SBJ list:"
read list_name

declare -a SBJs
SBJs=(`cat "${root_dir}scripts/SBJ_lists/${list_name}.sbj"`)

# Run BHV scripts for these SBJs
for sbj in "${SBJs[@]}"; do
    echo "=================================================================\n"
    echo "Running ${sbj}\n"
    echo "=================================================================\n"
    python BHV00_extract.py ${sbj}
#     if [ ! -f ${root_dir}/data/${sbj}/03_events/${sbj}_prdm_vars.pkl ]; then
#         echo "pkl of paradigm vars not found, running for ${sbj}..."
#         python BHV00_extract.py ${sbj}
#     else
#         echo "pkl found, continuing to BHV01 for ${sbj}..."
#     fi
    python BHV01_prelim_analysis.py ${sbj}
done

