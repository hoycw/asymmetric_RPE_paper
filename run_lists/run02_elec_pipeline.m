%% Check which root directory
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load SBJ list
SBJ_id = 'preproc';
SBJs = fn_load_SBJ_list(SBJ_id);

%% Prepare for inspections
% Requires SBJ_vars to be set up

for s = 1:numel(SBJs)
    % Convert raw pipeline elec files to my SBJ_vars
%     fn_elec_import_orig(SBJs{s},'main_ft','pat','',0);
%     fn_elec_import_orig(SBJs{s},'main_ft','pat','',1);
%     fn_elec_import_orig(SBJs{s},'main_ft','mni','v',1);
    
    % Match elec to atlas labels + tissue (ONLY orig!)
    % run in SGE: fn_elec_match_atlas(SBJs{s},'main_ft','pat','','Dx');
    
    % Export reref atlas info to CSV for manual adjustments
%     fn_elec_export_csv(SBJs{s},'main_ft','pat','','Dx', 0);
    
end

%% Manual adjustments of orig elec
% Import .csv to google sheet, then check assignments for manual adjustments
%   -duplicate to keep coloring
%   -import new data at cell 2A
%       -File -> Import; Replace data at selected cell; convert text to numbers, data, and formulas
%   -delete remaining rows from old data

% Run subsections of this manually to plot elecs relative to pial, white
%   matter, and inflated surfaces for each probe
% fn_elec_check_ROIs(SBJ);

%% Compile manual orig into auto bipolar
% Download manual adjustments into .tsv
for s = 1:numel(SBJs)
    fn_elec_import_manual(SBJs{s}, 'main_ft', 'pat', '', 'Dx', 0);
    fn_elec_compile_man_reref(SBJs{s}, 'main_ft', 'pat', '', 'Dx');
    fn_elec_export_csv(SBJs{s},'main_ft','pat','','Dx', 1);
end

%% Manual adjustments of reref elec
% Import .csv to google sheet, then check assignments for manual adjustments
%   -especially check logic for merging ROIs

% fn_elec_check_ROIs(SBJ);

%% Complete manual updates
for s = 1:numel(SBJs)
    % Reimport manual adjustments
    fn_elec_import_manual(SBJs{s}, 'main_ft', 'pat', '', 'Dx', 1);
    
    % Copy corrected labels to MNI elec files
    fn_elec_copy_atlas_pat2mni(SBJs{s},'main_ft','v','Dx');
end

%% MNI Check ???
% SBJ = 'IR21';
% 
% % Compare patient and MNI in ortho
% fn_view_recon(SBJ,'main_ft','ortho','pat','',1,'b',1);
% fn_view_recon(SBJ,'main_ft','ortho','mni','v',1,'b',1);
% 
% % Check atlas assignments
% fn_view_recon_atlas(SBJ,proc_id,'pat','',1,'b','DK','gROI');
% fn_view_recon_atlas(SBJ,proc_id,'pat','',1,'b','Dx','gROI');
% fn_view_recon_atlas(SBJ,proc_id,'mni','v',1,'b','Yeo7','Yeo7');

%% ================================================================================
%  PLOT ANATOMY
%  =================================================================================
%% Plot group recon with mgROI
% fn_view_recon_atlas_grp(SBJs,proc_id,'v',0,'l','Dx','mgROI',0);
% fn_view_recon_atlas_grp(SBJs,proc_id,'v',0,'r','Dx','mgROI',0);

roi_opts  = {{'l','INS',1},{'l','lat',1},{'l','MPFC',1}};%,{'b','OFC',0}};,{'l','deep',1}};
proc_id   = 'main_ft';
roi_id    = 'main3';%'mgROI';
atlas_id  = 'Dx';
reg_type  = 'v';
show_lab  = 0;
save_fig  = 1;
fig_ftype = 'png';

for roi_ix = 1:numel(roi_opts)
    fn_view_recon_atlas_grp_ROI(SBJ_id, proc_id, reg_type, show_lab,...
                                roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
                                roi_opts{roi_ix}{3},'save_fig', save_fig, 'fig_ftype', fig_ftype);
end

%% Plot Single Electrode
SBJ = 'IR57';
plot_elecs = {'RAC1-2'};
roi_opts  = {{'r','MPFC'}};
proc_id   = 'main_ft';
atlas_id  = 'Dx';

SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final.mat']);
for e = 1:numel(plot_elecs)
    % Select and color electrode
    cfgs = []; cfgs.channel = plot_elecs(e);
    ch_elec = fn_select_elec(cfgs,elec);
    ch_elec.color = fn_roi2color(ch_elec.gROI);
    view_angle = fn_get_view_angle(roi_opts{e}{1},roi_opts{e}{2});

    % Load Atlas
    atlas = fn_load_recon_atlas([],atlas_id);
    atlas_labels = fn_atlas_roi_select_mesh(atlas_id, roi_opts{e}{2}, roi_opts{e}{1});
    
    % Create ROI mesh
    cfg = [];
    cfg.inputcoord = atlas.coordsys;
    cfg.atlas = atlas;
    cfg.roi = atlas_labels;
    roi_mask = ft_volumelookup(cfg,atlas);
    seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
    seg.brain = roi_mask;
    cfg = [];
    cfg.method      = 'iso2mesh';   % surface toolbox Arjen found
    cfg.radbound    = 2;            % scalar indicating the radius of the target surface mesh element bounding sphere
    cfg.maxsurf     = 0;
    cfg.tissue      = 'brain';
    cfg.numvertices = 100000;
    cfg.smooth      = 3;
    cfg.spmversion  = 'spm12';
    roi_mesh = ft_prepare_mesh(cfg, seg);
    
    % Plot Recon
    figure;
    ft_plot_mesh(roi_mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.3);
    ft_plot_sens(ch_elec, 'elecshape', 'sphere', 'facecolor', ch_elec.color, 'label', 'off');
    view(view_angle); material dull; lighting gouraud;
    l = camlight;
end

%% Print MPFC electrodes
for s = 1:numel(SBJs)
%     SBJ = SBJs{s};
%     fprintf('%s\n',SBJ);
%     load([root_dir 'PRJ_Error/data/' SBJ '/05_recon/' SBJ '_elec_main_ft_pat_Dx_full.mat']);
%     if any(strcmp(SBJ,{'CP24','IR57','IR68'}))
%         elec.roi = fn_atlas2roi_labels(elec.atlas_lab,'Dx','gROI');
%     end
%     mpfc_ix = find(strcmp(elec.roi,'MPFC'));
%     disp(elec.label(mpfc_ix));
end

%% Manual vs. Automatic ROI comparison
roi_id = 'ROI';
atlas_id = 'Dx';
proc_id = 'main_ft';

% SBJs = {'CP24','IR21','IR26','IR31','IR32','IR35','IR39','IR41',...
%         'IR52','IR54','IR57','IR61','IR65','IR67','IR68','IR72','IR74'};
% 
% roi_list = fn_roi_label_styles(roi_id);
% change = zeros(numel(roi_list));
% roi_cnt = zeros([numel(roi_list) 1]);
% n_elec = 0;
% for s = 1:numel(SBJs)
%     eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJs{s} '_vars.m']);
%     fin_fname = [SBJ_vars.dirs.recon,SBJs{s},'_elec_',proc_id,'_pat','_',atlas_id,'_final.mat'];
%     man_fname = [SBJ_vars.dirs.recon,SBJs{s},'_elec_',proc_id,'_pat','_',atlas_id,'_man.mat'];
%     load(fin_fname); fin = elec;
%     load(man_fname); man = elec;
%     if numel(fin.label)~=numel(man.label); error('num chan off');end
%     n_elec = n_elec + numel(man.label);
%     for e = 1:numel(fin.label)
%         fin_ix = find(strcmp(roi_list,fin.(roi_id){e}));
%         man_ix = find(strcmp(roi_list,man.(roi_id){e}));
%         change(man_ix,fin_ix) = change(man_ix,fin_ix) + 1;
%     end
% end
% change_norm = zeros(numel(roi_list));
% change_z = change;
% for r = 1:numel(roi_list)
%     change_norm(r,:) = change(r,:)/sum(change(r,:));
%     change_z(r,r) = 0;
% end

%%
% figure;
% imagesc(change_norm);set(gca,'YDir','normal');
% colorbar;
% set(gca,'XTick',1:1:numel(roi_list));
% set(gca,'XTickLabel',roi_list);
% set(gca,'YTick',1:1:numel(roi_list));
% set(gca,'YTickLabel',roi_list);
% title('normalized');
% 
% figure;
% imagesc(change_z);set(gca,'YDir','normal');
% colorbar;
% set(gca,'XTick',1:1:numel(roi_list));
% set(gca,'XTickLabel',roi_list);
% set(gca,'YTick',1:1:numel(roi_list));
% set(gca,'YTickLabel',roi_list);
% title('count without diagonal');
