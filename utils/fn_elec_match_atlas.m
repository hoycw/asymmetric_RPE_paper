function fn_elec_match_atlas(SBJ, proc_id, view_space, reg_type, atlas_id)
%% Match elec_orig to atlas ROIs and tissue
% INPUTS:
%   SBJ [str] - subject ID to plot
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   view_space [str] - {'pat', 'mni'}
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'} are the only ones implemented so far

if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];
else
    reg_suffix = '';
end

%% Load elec struct
load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_',view_space,reg_suffix,'_orig.mat']);

%% Load Atlas
atlas = fn_load_recon_atlas(SBJ,atlas_id);

%% Match elecs to atlas ROIs
elec = fn_atlas_lookup(elec,atlas,'min_qry_rng',1,'max_qry_rng',5);

%% Match elecs to atlas tissue compartments
elec.roi_flag   = zeros(size(elec.label));
if any(strcmp(atlas_id,{'DK','Dx'}))
    tiss = fn_atlas_lookup(elec,atlas,'min_qry_rng',5,'max_qry_rng',5);
    
    %% Convert atlas labels and probabilities to GM probability
    % usedqueryrange search sizes: 1 = 1; 3 = 7; 5 = 33
    elec.tissue_labels = {'GM','WM','CSF','OUT'};
    elec.tissue_prob = zeros([numel(elec.label) numel(elec.tissue_labels)]);
    
    % Assign atlas labels to tissue type
    elec.tissue    = fn_atlas2roi_labels(tiss.atlas_lab,atlas_id,'tissue');
    elec.tissue2   = cell(size(elec.tissue));
    elec.gm_weight = zeros(size(elec.label));
    for e = 1:numel(elec.label)
        % Compute Probability of Tissue Types {GM, WM, CSF, OUT}
        elec.tissue_prob(e,strcmp(elec.tissue{e},elec.tissue_labels)) = tiss.atlas_prob(e);
        
        % Check for secondary matches and add to total
        if ~isempty(tiss.atlas_lab2{e})
            elec.tissue2{e} = fn_atlas2roi_labels(tiss.atlas_lab2{e},atlas_id,'tissue');
            for roi = 1:numel(elec.tissue2{e})
                elec.tissue_prob(e,strcmp(elec.tissue2{e}{roi},elec.tissue_labels)) = ...
                    elec.tissue_prob(e,strcmp(elec.tissue2{e}{roi},elec.tissue_labels)) + tiss.atlas_prob2{e}(roi);
            end
        end
        
        % Add gm_weight
        elec.gm_weight(e) = fn_gm_weight(elec.tissue{e},elec.tissue2{e});
        
        % Add roi_flag for non-GM/WM
        if any(strcmp(elec.tissue{e},{'CSF','OUT'}))
            elec.roi_flag(e) = 1;
        end
    end
end

%% Add ROI and gROI
elec.gROI = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,'gROI');
elec.ROI  = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,'ROI');

%% Save elec strcut with atlas labels
out_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_',view_space,reg_suffix,'_orig_',atlas_id,'.mat'];
fprintf('Saving %s\n',out_fname);
fprintf('==================================================================\n');
save(out_fname,'-v7.3','elec');

end
