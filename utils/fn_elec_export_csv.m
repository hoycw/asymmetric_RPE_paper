function fn_elec_export_csv(SBJ, proc_id, view_space, reg_type, atlas_id, reref)
%% Export elec struct to CSV for manual edits on google spreadsheet
%   If orig, adds gROI and ROI

%% Load elec and get ROI labels
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
else
    elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space reg_suffix '_orig_' atlas_id '.mat'];
end
load(elec_fname);

%% Print data to csv file
csv_fname = [elec_fname(1:end-4) '.csv'];
fprintf('\tExporting %s...\n',csv_fname);
csv = fopen(csv_fname,'w');

for e = 1:numel(elec.label)
    if ~reref
        if ~isempty(elec.atlas_lab2{e})
            if ischar(elec.atlas_lab2{e})
                atlas_lab2 = elec.atlas_lab2{e};
            else
                atlas_lab2 = strjoin(elec.atlas_lab2{e},'+');
            end
        else
            atlas_lab2 = '';
        end
        % Export order:
        %   label, atlas_lab, atlas_prob, atlas_lab2, atlas_qryrng, ...
        %   tissue, tissue_prob (4x), gm_weight, hemi, gROI, ROI, roi_flag
        fprintf(csv,'%s,%s,%.3f,%s,%d,%s,%.3f,%.3f,%.3f,%.3f,%.1f,%s,%s,%s,%d,0,0,0,%s\n',...
            elec.label{e},elec.atlas_lab{e},elec.atlas_prob(e),atlas_lab2,elec.atlas_qryrng(e),...
            elec.tissue{e},elec.tissue_prob(e,:),elec.gm_weight(e),elec.hemi{e},...
            elec.gROI{e},elec.ROI{e},elec.roi_flag(e),'');
    else
        % Export order:
        %   label, ROI1, gm_weight1, par_vol1, ROI2, gm_weight2, par_vol2,
        %   tissue, man_adj, hemi, gROI, ROI, roi_flag, re_adj (0), border (0), anat_notes
        fprintf(csv,'%s,%s,%s,%d,%s,%s,%d,%s,%.1f,%s,%s,%s,%d,0,0,%s\n',...
            elec.label{e},elec.inputs{e}.ROI{1},elec.inputs{e}.tissue{1},elec.inputs{e}.par_vol(1),...
            elec.inputs{e}.ROI{2},elec.inputs{e}.tissue{2},elec.inputs{e}.par_vol(2),elec.tissue{e},...
            elec.man_adj(e),elec.hemi{e},elec.gROI{e},elec.ROI{e},elec.roi_flag(e),elec.anat_notes{e});
    end
end

fclose(csv);

end

