%% Run scripts for SfN 2019 abstract
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath(genpath([root_dir 'PRJ_Error/scripts/']));
addpath(genpath([root_dir 'PRJ_Error/scripts/utils/']));
addpath(ft_dir);
ft_defaults

%%
behav_SBJ = {'CP22','CP23','CP24','CP25','IR57','IR60','IR62','IR63','IR65','IR66','IR67','IR68',...
                'IR69','IR71','IRR72','IR74','IR75','IR76','IR77','IR78','IR79','IR82','IR84'};
eeg_behav_SBJ = {'EP06','EP07','EP08','EP10','EP11','EP14','EP15','EP16','EP17','EP18','EP19',...
           'EEG01','EEG02','EEG03','EEG04','EEG06','EEG07','EEG08','EEG09','EEG10','EEG12'};
SBJs = {'CP24','IR57','IR68'};
proc_id    = 'main_ft';
actv_win   = 100;
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

%% PLOT GROUP ROI MESHES
roi_opts  = {{'l','deep',1},{'l','lat',1},{'l','MPFC',1},{'b','OFC',0}};
%,{'r','deep'},{'r','lat'},{'r','MPFC'}
proc_id   = 'main_ft';
atlas_id  = 'Dx';
roi_id    = 'mgROI';
reg_type  = 'v';
show_lab  = 0;
save_fig  = 1;
fig_ftype = 'png';

for roi_ix = 2:numel(roi_opts)
    fn_view_recon_atlas_grp_ROI(SBJs, proc_id, reg_type, show_lab,...
                                roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
                                roi_opts{roi_ix}{3},'save_fig', save_fig, 'fig_ftype', fig_ftype);
end

%% Plot Group Stat ROI Recon
proc_id   = 'main_ft';
an_id   = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';
stat_id = 'corrRT_DifOutTimDO_WL200_WS50';%'DifOutDO_F0t1_WL1_WS25';
hemi        = 'l';
roi_id      = 'mgROI';
atlas_id    = 'Dx';
reg_type    = 'v';
show_labels = 0;
save_fig    = 1;
fig_ftype   = 'png';

roi_opts  = {{'l','deep',1},{'l','lat',1},{'l','MPFC',1},{'b','OFC',0}};
for roi_ix = 1:numel(roi_opts)
    fn_view_recon_atlas_grp_stat_ROI(SBJs, proc_id, stat_id, an_id, ...
        reg_type, show_labels, roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
        roi_opts{roi_ix}{3},'save_fig', save_fig, 'fig_ftype', fig_ftype);
    fn_view_recon_atlas_grp_stat_venn_ROI(SBJs, proc_id, stat_id, an_id,...
        reg_type, show_labels, roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
        roi_opts{roi_ix}{3},'save_fig', save_fig, 'fig_ftype', fig_ftype);
end

%%
an_id   = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';
stat_id = 'DifOutDO_F0t1_WL1_WS25';
for s = 1:numel(SBJs)
%     SBJ08a_crRT_nANOVA(SBJs{s},an_id,stat_id);
end

%%
% conditions = 'DifOut';
proc_id   = 'main_ft';
an_id     = 'HGm_S_zbtS_trl2to3001_sm0_wn100_stat3';%'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';%
stat_id   = 'corrRT_DifOutTimDO_WL200_WS50';
plt_id    = 'ts_S0to3_evnts_sigline';%'ts_F0to1_evnts_sigline';
atlas_id  = 'Dx';
roi_id    = 'mgROI';
gm_thresh = 0;
z_thresh  = 0;
plot_nsig = 0;
save_fig  = 1;
fig_vis   = 'on';
for s = 1:numel(SBJs)
%     plt_id     = 'stack_S2to3_evnt_c5';%'stack_F2to12_evnt_c5';
%     SBJ07b_HFA_plot_stack_cond_saved(SBJs{s}, conditions, an_id, actv_win, plt_id, save_fig, fig_vis, fig_ftype);
%     close all;
    
%     SBJ08b_HFA_plot_corrRT_ANOVA(SBJs{s}, proc_id, an_id, stat_id, plt_id, save_fig, fig_vis, fig_ftype);
%     close all;
    SBJ10b_ANOVA_plot_ts_gROIcomb(SBJs{s},proc_id,stat_id,an_id,...
                                    atlas_id,roi_id,gm_thresh,z_thresh,plot_nsig,...
                                    plt_id,save_fig,fig_vis,fig_ftype)
end

%%
roi_id = 'gROI';
atlas_id = 'Dx';
gm_thresh = 0;
plot_out = 0;
plot_scat = 1;
plt_id      = '';
fig_ftype = 'png';
SBJ08c_HFA_GRP_summary_errbar_perc_GMlim_actv_RT_ANOVA_ROI(SBJs,stat_id,proc_id,an_id,actv_win,roi_id,...
                                                            atlas_id,gm_thresh,plt_id,plot_out,plot_scat,save_fig,fig_vis,fig_ftype);
SBJ08c_HFA_GRP_summary_errbar_perc_GMlim_actv_RT_ANOVA_ROI_DO(SBJs,stat_id,proc_id,an_id,actv_win,roi_id,...
                                                            atlas_id,gm_thresh,plt_id,plot_out,plot_scat,save_fig,fig_vis,fig_ftype);

%%
for s = 1:numel(SBJs)
    SBJ = SBJs{s};
    fprintf('%s\n',SBJ);
    load([root_dir 'PRJ_Error/data/' SBJ '/05_recon/' SBJ '_elec_main_ft_pat_Dx_full.mat']);
    if any(strcmp(SBJ,{'CP24','IR57','IR68'}))
        elec.roi = fn_atlas2roi_labels(elec.atlas_lab,'Dx','gROI');
    end
    mpfc_ix = find(strcmp(elec.roi,'MPFC'));
    disp(elec.label(mpfc_ix));
end

%% Recon viewing
SBJs = {'CP24','IR57','IR68'};%'IR66' - doesn't have mni_v_Dx yet
an_id = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';%'HGm_S_zbtS_trl2to3001_sm0_wn100_stat3';%
mirror = 1;
for s = 1:numel(SBJs)
%     fn_view_recon_stat(SBJs{s},proc_id,stat_id,an_id, 'pat', '', 0, 'r', 0, 'save_fig', 1);
    fn_view_recon_stat(SBJs{s},proc_id,stat_id,an_id, 'pat', '', 0, 'l', mirror, 0, 'save_fig', 1);
    %close all;
end

%%
SBJs = {'CP24','IR57','IR68'};%'IR66' - doesn't have mni_v_Dx yet
roi_id = 'gROI';
atlas_id = 'Dx';
plot_out = 0;
show_lab = 0;
fn_view_recon_atlas_grp_stat(SBJs, proc_id, stat_id, an_id, 'v', show_lab, 'r', atlas_id, roi_id, plot_out);
fn_view_recon_atlas_grp_stat(SBJs, proc_id, stat_id, an_id, 'v', show_lab, 'l', atlas_id, roi_id, plot_out);