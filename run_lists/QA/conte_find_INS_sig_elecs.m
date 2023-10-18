%% Find specific ROI of RPE effects from thesis for Conte center update
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath([root_dir 'PRJ_Error/scripts/']);
addpath(genpath([root_dir 'PRJ_Error/scripts/utils/']));
addpath(ft_dir);
ft_defaults

%%
SBJ_id = 'preproc';%'preproc_nu';%
SBJs = fn_load_SBJ_list(SBJ_id);
proc_id   = 'main_ft';
an_id     = 'HGm_F25t121_zbtS_sm0_l1_wn50';
model_id  = 'ERPEs_DifFB';
stat_id   = 'mGLM_st0t6_WL05_WS25';
atlas_id  = 'Dx';
roi_id    = 'INS';
hemi      = 'b';

eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);
[reg_lab, reg_names, reg_colors, ~, ~] = fn_regressor_label_styles(mdl.model_lab);

%%
elec_sbj    = cell([numel(SBJs) 1]);
ins_elec_cnt = zeros(size(SBJs));
sRPE_elecs = {};
sRPE_ROI   = {};
uRPE_elecs = {};
uRPE_ROI   = {};
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load elec
    elec_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final.mat'];
    tmp = load(elec_fname); elec_sbj{sbj_ix} = tmp.elec;
    ins_elec_cnt(sbj_ix) = sum(strcmp(elec_sbj{sbj_ix}.gROI,'INS'));
    
    % Append SBJ name to labels
    orig_labels = elec_sbj{sbj_ix}.label;
    for e_ix = 1:numel(elec_sbj{sbj_ix}.label)
        elec_sbj{sbj_ix}.label{e_ix} = [SBJs{sbj_ix} '_' elec_sbj{sbj_ix}.label{e_ix}];
    end
    roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, hemi, atlas_id, roi_id);
    
    %% Find significant electrodes in INS
    if ~isempty(roi_elecs)
        load([SBJ_vars.dirs.stats SBJ '_mGLM_ROI_' model_id '_' stat_id '_' an_id '.mat'],'beta');
        for e_ix = 1:numel(beta.label)
            sbj_e_lab = [SBJs{sbj_ix} '_' beta.label{e_ix}];
            % Add if sRPE is significant
            if any(squeeze(beta.qval(strcmp(beta.feature,'sRPE'),e_ix,:)<=0.05)) && ...
                    contains(sbj_e_lab,roi_elecs)
                sRPE_elecs = [sRPE_elecs; sbj_e_lab];
                
                sRPE_ROI = [sRPE_ROI; elec_sbj{sbj_ix}.ROI(strcmp(elec_sbj{sbj_ix}.label,sbj_e_lab))];
            end
            % Add if uRPE is significant
            if any(squeeze(beta.qval(strcmp(beta.feature,'uRPE'),e_ix,:)<=0.05)) && ...
                    contains(sbj_e_lab,roi_elecs)
                uRPE_elecs = [uRPE_elecs; sbj_e_lab];
                uRPE_ROI = [uRPE_ROI; elec_sbj{sbj_ix}.ROI(strcmp(elec_sbj{sbj_ix}.label,sbj_e_lab))];
            end
        end
    end
    
    clear SBJ SBJ_vars SBJ_vars_cmd beta
end

%%
% fn_view_recon('IR87','main_ft','ortho','pat','',1,'b',1);