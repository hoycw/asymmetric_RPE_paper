%% Run scripts for SfN 2019 abstract
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath(genpath([root_dir 'PRJ_Error/scripts/']));
addpath(genpath([root_dir 'PRJ_Error/scripts/utils/']));
addpath(ft_dir);
ft_defaults

%%
SBJs = {'CP24','IR57','IR68'};

%% Run Model (SGE...)
an_id   = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';%'HGm_F25t121_zbtS_sm0_l1_wn100';
stat_id = 'RL_DifOut_F01t06_WL005_WS005';
for s = 1:numel(SBJs)
    SBJ08a_crRT_mGLM(SBJs{s},an_id,stat_id);
end

%% Plot Example HFA stack and ANOVA ts
SBJ        = 'IR57';
elec_lab   = {'RAC1-2'};
an_id      = 'HGm_S_zbtS_trl2to3001_sm0_wn100_stat3';
conditions = 'DifOut';
actv_win   = 100;
plt_id     = 'stack_S2to3_evnt_c5';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'svg';
SBJ07b_HFA_plot_stack_cond_saved_single(SBJ, elec_lab, conditions, an_id, actv_win,...
    plt_id, save_fig, fig_vis, fig_ftype);

an_id      = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';
stat_id    = 'DifOutDO_F0t1_WL1_WS25';
plt_id    = 'ts_F0to1_evnts_sigline';
SBJ08b_HFA_plot_crRT_nANOVA_single(SBJs{s}, elec_lab, an_id, stat_id, plt_id, save_fig, fig_vis, fig_ftype);
close all;

%% Plot Group Stat ROI Recon
% proc_id   = 'main_ft';
% an_id     = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';
% stat_id   = 'DifOutDO_F0t1_WL1_WS25';%'corrRT_DifOutTimDO_WL200_WS50';%
% hemi      = 'l';
% roi_id    = 'mgROI';
% atlas_id  = 'Dx';
% reg_type  = 'v';
% show_lab  = 0;
% save_fig  = 1;
% fig_ftype = 'png';
% 
% roi_opts  = {{'l','deep',1},{'l','lat',1},{'l','MPFC',1},{'b','OFC',0}};
% for roi_ix = 1:numel(roi_opts)
% %     fn_view_recon_atlas_grp_stat_ROI(SBJs, proc_id, stat_id, an_id, ...
% %         reg_type, show_labels, roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
% %         roi_opts{roi_ix}{3},'save_fig', save_fig, 'fig_ftype', fig_ftype);
%     fn_view_recon_atlas_grp_stat_venn_ROI(SBJs, proc_id, stat_id, an_id,...
%         reg_type, show_lab, roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
%         roi_opts{roi_ix}{3},'save_fig', save_fig, 'fig_ftype', fig_ftype);
% end

%% Plot ANOVA time series
% % conditions = 'DifOut';
% proc_id   = 'main_ft';
% an_id     = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';%'HGm_S_zbtS_trl2to3001_sm0_wn100_stat3';%
% stat_id   = 'DifOutDO_F0t1_WL1_WS25';%'corrRT_DifOutTimDO_WL200_WS50';
% plt_id    = 'ts_F0to1_evnts_sigline';
% atlas_id  = 'Dx';
% roi_id    = 'mgROI';
% save_fig  = 1;
% fig_vis   = 'on';
% fig_ftype  = 'svg';
% SBJ10c_grp_ANOVA_plot_ts_gROIcomb(SBJs,proc_id,stat_id,an_id,...
%                                     atlas_id,roi_id,plt_id,save_fig,fig_vis,fig_ftype);

%%
% an_id     = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';
% stat_id   = 'DifOutDO_F0t1_WL1_WS25';
% atlas_id  = 'Dx';
% roi_id    = 'mgROI';
% gm_thresh = 0;
% plot_out  = 0;
% plot_scat = 1;
% plt_id    = '';
% save_fig  = 1;
% fig_vis   = 'on';
% fig_ftype = 'svg';
% SBJ08c_HFA_grp_errbar_ROI_nANOVA(SBJs,stat_id,proc_id,an_id,atlas_id,roi_id,...
%                                           plot_scat,save_fig,fig_vis,fig_ftype);
                                      