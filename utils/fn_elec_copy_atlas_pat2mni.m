function fn_elec_copy_atlas_pat2mni(SBJ,proc_id,reg_type,atlas_id)%,varargin)
%% Copy atlas+tissue info from patient space elec to mni space elec
% INPUTS:
%   SBJ [str] - name of subject
%   proc_id [str] - name of analysis pipeline
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}

% Set up paths
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

% %% Variable input handling
% if ~isempty(varargin)
%     for v = 1:2:numel(varargin)
%         if strcmp(varargin{v},'roi_id')
%             roi_id = varargin{v+1};
%         else
%             error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
%         end
%     end
% end

%% Load variables
eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/proc_vars/' proc_id '_vars.m']);

%% Load Elec struct
if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];
else
    error('reg_type must be selected for mni space');
end

pat_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat','_',atlas_id,'_final.mat'];
mni_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_mni',reg_suffix,'.mat'];
out_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_mni',reg_suffix,'_',atlas_id,'_final.mat'];
% if exist('roi_id','var')
%     if exist([pat_fname(1:end-4) '_' roi_id '.mat'],'file')
%         pat_fname = [pat_fname(1:end-4) '_' roi_id '.mat'];
%         out_fname = [out_fname(1:end-4) '_' roi_id '.mat'];
%     end
% end
load(pat_fname); pat = elec;
load(mni_fname);

% Check that all elecs are in both structs
if ~all(strcmp(elec.label,pat.label))
    elec = fn_reorder_elec(elec,pat.label);
    if ~all(strcmp(elec.label,pat.label))
        fprintf(2,'Missing channels in %s:\n',mni_fname);
        disp(setdiff(elec.label,pat.label));
        error('Fix these missing channels before copying atlas info');
    end
end

%% Add Back Stripped Fields Channel Types
fields = fieldnames(pat);
fields = setdiff(fields,fieldnames(elec));
% Overwrite fields that were manually editted
fields = [fields; {'tissue'; 'hemi'; 'gROI'; 'ROI'}];
% fields = {'atlas_lab', 'atlas_lab2', 'gm_weight', 'hemi', 'inputs',...
%             'roi_flag', 'tissue', 'tissue2', 'tissue_prob'};
for f = 1:numel(fields)
    elec.(fields{f}) = pat.(fields{f});
end

%% Save data
% Check if elec.cfg.previous got ridiculously large, and keep only first
var_stats = whos('elec');
if var_stats.bytes>1000000
    elec.cfg = rmfield(elec.cfg,'previous');
end
fprintf('============== Saving %s ==============\n',out_fname);
save(out_fname, '-v7.3', 'elec');

end
