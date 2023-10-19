#!/bin/sh

# Move to scripts directory 
root_dir="/Volumes/hoycw_clust/PRJ_Error"
out_dir="/Users/colinhoy/Code/PRJ_Error/asymmetric_RPE_paper_data"
cd $root_dir/scripts/

# Load the correct SBJ list
echo "Please enter the desired SBJ list:"
read list_name

declare -a SBJs
SBJs=(`cat "${root_dir}/scripts/SBJ_lists/${list_name}.sbj"`)

# Run BHV scripts for these SBJs
for sbj in "${SBJs[@]}"; do
    echo "=================================================================\n"
    echo "Running ${sbj}\n"
    echo "=================================================================\n"
    mkdir ${out_dir}/${sbj}/
    mkdir ${out_dir}/${sbj}/02_preproc/
    mkdir ${out_dir}/${sbj}/03_events/
    mkdir ${out_dir}/${sbj}/04_proc/
    mkdir ${out_dir}/${sbj}/05_recon/
    cp ${root_dir}/data/${sbj}/02_preproc/${sbj}_preproc_main_ft.mat ${out_dir}/${sbj}/02_preproc/
    cp ${root_dir}/data/${sbj}/03_events/${sbj}_bhv_main_ft_final.mat ${out_dir}/${sbj}/03_events/
    cp ${root_dir}/data/${sbj}/03_events/${sbj}_prdm_vars.mat ${out_dir}/${sbj}/03_events/
    cp ${root_dir}/data/${sbj}/04_proc/${sbj}_ROI_main_ft_HGm_F25t121_zbtS_sm0_l1_wn50.mat ${out_dir}/${sbj}/04_proc/
    cp ${root_dir}/data/${sbj}/05_recon/${sbj}_elec_main_ft_mni_v_Dx_final.mat ${out_dir}/${sbj}/05_recon/
done

