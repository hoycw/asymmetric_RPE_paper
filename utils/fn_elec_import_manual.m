function fn_elec_import_manual(SBJ, proc_id, view_space, reg_type, atlas_id, reref)
%% Import TSV file with manual adjustments of elec files
%   Added field:
%   elec.man_adj [0/1] - was this adjusted manually
%   elec.notes [str] - alternative ROIs usually

%% Load SBJ parameters and elec
% Check which root directory
[root_dir, ~] = fn_get_root_dir();
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];
else
    reg_suffix = '';
end

if reref
    elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space reg_suffix '_' atlas_id '_man.mat'];
    tsv_fname  = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space reg_suffix '_' atlas_id '_final.tsv'];
    save_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space reg_suffix '_' atlas_id '_final.mat'];
else
    elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space reg_suffix '_orig_' atlas_id '.mat'];
    tsv_fname  = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space reg_suffix '_orig_' atlas_id '_man.tsv'];
    save_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space reg_suffix '_orig_' atlas_id '_man.mat'];
end
load(elec_fname);

% %% Add ROI fields for orig elecs
% if ~reref
%     elec.gROI = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,'gROI');
%     elec.ROI  = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,'ROI');
% end
% 
%% Load manually adjusted TSV
tsv_file = fopen(tsv_fname, 'r');
if ~reref
    % Original Elecs
    headers = textscan(tsv_file, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 1);
    % Import order:
    %   label, atlas_lab, atlas_prob, atlas_lab2, atlas_qryrng, ...
    %   tissue, tissue_prob (4x), gm_weight, hemi, gROI, ROI, roi_flag,
    %   man_adj, par_vol, WM_ref, anatomy notes
    man_adj = textscan(tsv_file, '%s %s %f %s %d %s %f %f %f %f %f %s %s %s %d %d %d %d %s', 'HeaderLines', 1,...
        'Delimiter', '\t', 'MultipleDelimsAsOne', 0);
    new_fields = {'tissue','hemi','gROI','ROI','man_adj','par_vol','WM_ref','anat_notes'};
else
    % Rereferenced Elecs
    headers = textscan(tsv_file, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 1);
    % Import order:
    %   label, ROI1, gm_weight1, par_vol1, ROI2, gm_weight2, par_vol2, ...
    %   tissue, man_adj, hemi, gROI, ROI, roi_flag, re_adj, anat_notes
    man_adj = textscan(tsv_file, '%s %s %s %f %s %s %f %s %f %s %s %s %d %d %d %s', 'HeaderLines', 1,...
        'Delimiter', '\t', 'MultipleDelimsAsOne', 0);
    new_fields = {'tissue','hemi','gROI','ROI','re_adj','border','anat_notes'};
end
fclose(tsv_file);

%% Add data back to elec
for f = 1:numel(new_fields)
    man_adj_ix = find(strcmp([headers{:}], new_fields{f}));
    if strcmp(new_fields{f},'man_adj')
        elec.man_adj = man_adj{man_adj_ix};
    elseif strcmp(new_fields{f},'re_adj')
        elec.re_adj = man_adj{man_adj_ix};
    elseif strcmp(new_fields{f},'border')
        elec.border = man_adj{man_adj_ix};
    elseif strcmp(new_fields{f},'anat_notes')
        elec.anat_notes = man_adj{man_adj_ix};
        % If anat_notes is empty, last line isn't read
        if numel(elec.anat_notes)==numel(elec.label)-1
            elec.anat_notes{numel(elec.label)} = '';
        elseif numel(elec.anat_notes)~=numel(elec.label)
            error('anat_notes isnt matching up');
        end
    else
        elec.(new_fields{f}) = man_adj{man_adj_ix};
    end
end

%% Save new elec
fprintf('Saving %s\n',save_fname);
save(save_fname,'-v7.3','elec');

end

