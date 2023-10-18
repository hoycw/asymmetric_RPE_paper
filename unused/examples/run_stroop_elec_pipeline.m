%% Check which root directory
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Full SBJ list (completed and ready to run)
SBJs = {'IR68'};

%% Prepare for inspections
for s = 1:numel(SBJs)
    % Convert raw pipeline elec files to my SBJ_vars
%     fn_elec_import_orig(SBJs{s},'main_ft','pat','',0);
%     fn_elec_import_orig(SBJs{s},'main_ft','pat','',1);
%     fn_elec_import_orig(SBJs{s},'main_ft','mni','v',1);
    
    % Match elec to atlas labels + tissue (ONLY orig!)
    % run in SGE: fn_elec_match_atlas(SBJs{s},'main_ft','pat','','Dx');
    
    % Export reref atlas info to CSV for manual adjustments
    fn_elec_export_csv(SBJs{s},'main_ft','pat','','Dx', 0);
    
end

%% Manual adjustments of orig elec
% fn_elec_check_ROIs(SBJ);

%% Compile manual orig into auto bipolar
for s = 1:numel(SBJs)
    fn_elec_import_manual(SBJs{s}, 'main_ft', 'pat', '', 'Dx', 0);
    fn_elec_compile_man_reref(SBJs{s}, 'main_ft', 'pat', '', 'Dx');
    fn_elec_export_csv(SBJs{s},'main_ft','pat','','Dx', 1);
end

%% Complete manual updates
for s = 1:numel(SBJs)
    % Reimport manual adjustments
%     fn_elec_import_manual(SBJs{s}, 'main_ft', 'pat', '', 'Dx', 1);
    
    % Copy corrected labels to MNI elec files
    fn_elec_copy_atlas_pat2mni(SBJs{s},'main_ft','v','Dx');
end

%% MNI Check ???
% SBJ = 'IR21';
% 
% % Compare patient and MNI in ortho
% fn_view_recon(SBJ,'main_ft','ortho','pat','',1,'b');
% fn_view_recon(SBJ,'main_ft','ortho','mni','v',1,'b');
% 
% % Check atlas assignments
% fn_view_recon_atlas(SBJ,pipeline_id,'pat','',1,'b','DK','gROI');
% fn_view_recon_atlas(SBJ,pipeline_id,'pat','',1,'b','Dx','gROI');
% fn_view_recon_atlas(SBJ,pipeline_id,'mni','v',1,'b','Yeo7','Yeo7');

%% ROI comparison
roi_id = 'ROI';
atlas_id = 'Dx';
proc_id = 'main_ft';

SBJs = {'CP24','IR21','IR26','IR31','IR32','IR35','IR39','IR41',...
        'IR52','IR54','IR57','IR61','IR65','IR67','IR68','IR72','IR74'};

roi_list = fn_roi_label_styles(roi_id);
change = zeros(numel(roi_list));
roi_cnt = zeros([numel(roi_list) 1]);
n_elec = 0;
for s = 1:numel(SBJs)
    eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJs{s} '_vars.m']);
    fin_fname = [SBJ_vars.dirs.recon,SBJs{s},'_elec_',proc_id,'_pat','_',atlas_id,'_final.mat'];
    man_fname = [SBJ_vars.dirs.recon,SBJs{s},'_elec_',proc_id,'_pat','_',atlas_id,'_man.mat'];
    load(fin_fname); fin = elec;
    load(man_fname); man = elec;
    if numel(fin.label)~=numel(man.label); error('num chan off');end
    n_elec = n_elec + numel(man.label);
    for e = 1:numel(fin.label)
        fin_ix = find(strcmp(roi_list,fin.(roi_id){e}));
        man_ix = find(strcmp(roi_list,man.(roi_id){e}));
        change(man_ix,fin_ix) = change(man_ix,fin_ix) + 1;
    end
end
change_norm = zeros(numel(roi_list));
change_z = change;
for r = 1:numel(roi_list)
    change_norm(r,:) = change(r,:)/sum(change(r,:));
    change_z(r,r) = 0;
end

%%
figure;
imagesc(change_norm);set(gca,'YDir','normal');
colorbar;
set(gca,'XTick',1:1:numel(roi_list));
set(gca,'XTickLabel',roi_list);
set(gca,'YTick',1:1:numel(roi_list));
set(gca,'YTickLabel',roi_list);
title('normalized');

figure;
imagesc(change_z);set(gca,'YDir','normal');
colorbar;
set(gca,'XTick',1:1:numel(roi_list));
set(gca,'XTickLabel',roi_list);
set(gca,'YTick',1:1:numel(roi_list));
set(gca,'YTickLabel',roi_list);
title('count without diagonal');
