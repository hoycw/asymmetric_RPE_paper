%% Run scripts to get data for Thesis Table 1
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath([root_dir 'PRJ_Error/scripts/']);
addpath(genpath([root_dir 'PRJ_Error/scripts/utils/']));
addpath(ft_dir);
ft_defaults

%%
SBJ_id = 'preproc';
SBJs = fn_load_SBJ_list(SBJ_id);

proc_id = 'main_ft';
roi_id  = 'main3';
atlas_id = 'Dx';
[roi_list, ~, roi_field] = fn_roi_label_styles(roi_id);
[cond_lab, cond_names, cond_colors, ~, ~] = fn_condition_label_styles('DifFB');

%% Get Behavioral and Electrode Information
roi_cnt = zeros([numel(roi_list) numel(SBJs)]);
for s = 1:numel(SBJs)
    SBJ = SBJs{s};
    fprintf('%s:\n',SBJ);
    eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJs{s} '_vars.m']);
    load([root_dir 'PRJ_Error/data/' SBJ '/03_events/' SBJ '_bhv_' proc_id '_final.mat']);
    ez_idx = fn_condition_index({'Ez'}, bhv);
    s_idx = fn_condition_index({'Nu'},bhv);
    training_idx = bhv.blk==0;
    
    % Load full behavioral info
%     bhv_orig_tmp = cell(size(SBJ_vars.block_name));
%     for b_ix = 1:numel(SBJ_vars.block_name)
%         [bhv_orig_tmp{b_ix}] = fn_load_behav_csv([root_dir 'PRJ_Error/data/' SBJ '/03_events/' SBJ '_behav' SBJ_vars.block_name{b_ix} '.csv']);
%     end
%     
%     ez_idx_orig = fn_condition_index({'Ez'}, bhv_orig);
    
    % ITIs
    fprintf('\tITIs: '); fprintf('%.2f, ',unique(bhv.ITI_type)); fprintf('\n');
    
    % Number of trials
    fprintf('\tn_ez = %d; n_hd = %d\n',sum(ez_idx),sum(~ez_idx));
    
    % Accuracy (exclude surprise and training trials)
%     s_idx_orig = fn_condition_index({'Su'},bhv_orig);
    
    fprintf('\tez_acc = %.3f\n',mean(bhv.hit(ez_idx & ~s_idx & ~training_idx)));
    fprintf('\thd_acc = %.3f\n',mean(bhv.hit(~ez_idx & ~s_idx & ~training_idx)));

    % Electrodes
    elec_fname = [root_dir 'PRJ_Error/data/' SBJ '/05_recon/' SBJ '_elec_' proc_id '_pat_' atlas_id '_final.mat'];
    load(elec_fname);
    for roi_ix = 1:numel(roi_list)
        roi_cnt(roi_ix,s) = sum(strcmp(elec.(roi_field),roi_list{roi_ix}));
        fprintf('\t%s count = %d\n',roi_list{roi_ix},roi_cnt(roi_ix,s));
    end
    
    clear SBJ SBJ_vars bhv
end
