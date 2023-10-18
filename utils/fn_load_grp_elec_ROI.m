function [elec_sbj, good_sbj] = fn_load_grp_elec_ROI(SBJs,proc_id,atlas_id,roi_id,rcn)
%% Load elec files for each SBJ and prepare for mesh recon plotting
%   For each SBJ, load MNI space elec, append SBJ to labels, select elecs
%   within plot_roi and hemisphere, and mirror if necessary
% INPUTS:
%   SBJs [cell array] - list of SBJs to process
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   atlas_id [str] - ID of atlas to select ROIs: {'DK','Dx','Yeo7','Yeo17'}
%   roi_id [str] - ROI grouping by which to color the atlas ROIs
%   rcn [struct] - options for recon plotting, see fn_process_recon_vars.m
%       must contain .plot_roi_list, .mirror, .hemi, .reg_suffix
% OUTPUTS:
%   elec_sbj [Nx1 cell array] - elec struct for each SBJ
%   good_sbj [Nx1 bool array] - true for each SBJ excep tif no elecs remain

[root_dir, ~] = fn_get_root_dir();
[~, ~, roi_field] = fn_roi_label_styles(roi_id);

elec_sbj    = cell([numel(SBJs) 1]);
good_sbj = true(size(SBJs));
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('================= Processing elec for %s =================\n',SBJ);
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load elec
    elec_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_mni',rcn.reg_suffix,'_',atlas_id,'_final.mat'];
    tmp = load(elec_fname); elec_sbj{sbj_ix} = tmp.elec;

    % Append SBJ name to labels
    for e_ix = 1:numel(elec_sbj{sbj_ix}.label)
        elec_sbj{sbj_ix}.label{e_ix} = [SBJ ' ' elec_sbj{sbj_ix}.label{e_ix}];
    end
    
    % Match ROIs to colors
    elec_sbj{sbj_ix}.color = fn_roi2color(elec_sbj{sbj_ix}.(roi_field));

    % Find electrodes within atlas ROIs & hemisphere
    if rcn.mirror
        roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, 'b', atlas_id, roi_id);
    else
        roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, rcn.hemi, atlas_id, roi_id);
    end

    % Find elecrodes specifically within ROIs on plot_roi_list for mesh plotting
    plot_elecs = zeros([numel(elec_sbj{sbj_ix}.label) numel(rcn.plot_roi_list)]);
    for roi_ix = 1:numel(rcn.plot_roi_list)
        plot_elecs(:,roi_ix) = strcmp(elec_sbj{sbj_ix}.(roi_field),rcn.plot_roi_list{roi_ix});
    end

    % Select electrodes only within both categories
    cfgs = [];
    cfgs.channel = intersect(roi_elecs, elec_sbj{sbj_ix}.label(any(plot_elecs,2)));
    elec_sbj{sbj_ix} = fn_select_elec(cfgs, elec_sbj{sbj_ix});

    % Report outcome (and flip hemispheres if mirroring)
    if ~isempty(elec_sbj{sbj_ix}.label)
        fprintf('\t%s has %i channels in %s (hemi %s)\n',SBJ,size(elec_sbj{sbj_ix}.label,1),rcn.plot_roi,rcn.hemi_str);
        % Mirror hemispheres
        if rcn.mirror
            hemi_match_idx = strcmp(elec_sbj{sbj_ix}.hemi,rcn.hemi);
            if strcmp(rcn.plot_roi,'MPFCINS')
                ins_idx = strcmp(elec_sbj{sbj_ix}.gROI,'INS');
                % Flip INS to match rcn.hemi
                elec_sbj{sbj_ix}.chanpos(ins_idx & ~hemi_match_idx,1) = ...
                    -elec_sbj{sbj_ix}.chanpos(ins_idx & ~hemi_match_idx,1);
                % Flip MPFC to match opposite of rcn.hemi
                elec_sbj{sbj_ix}.chanpos(~ins_idx & hemi_match_idx,1) = ...
                    -elec_sbj{sbj_ix}.chanpos(~ins_idx & hemi_match_idx,1);
            else
                elec_sbj{sbj_ix}.chanpos(~hemi_match_idx,1) = ...
                    -elec_sbj{sbj_ix}.chanpos(~hemi_match_idx,1);
            end
        end
    else
        good_sbj(sbj_ix) = false;
        % Print no ROI match
        fprintf(2,'\t%s has no channels in %s ROI (hemi %s)\n',SBJ,rcn.plot_roi,rcn.hemi_str);
    end
    
    clear SBJ SBJ_vars SBJ_vars_cmd tmp plot_elecs roi_elecs cfgs
end

end