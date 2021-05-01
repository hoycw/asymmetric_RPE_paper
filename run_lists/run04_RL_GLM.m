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
stat_id = 'ERPEs_DO_glm_F0t5_WL01_WS25';
for s = 1:numel(SBJs)
    SBJ08a_crRT_mGLM(SBJs{s},an_id,stat_id);
end

%% Plot GLM ts
SBJ        = 'IR68';
proc_id    = 'main_ft';
an_id      = 'HGm_S_zbtS_trl2to3001_sm0_wn100_stat3';
conditions = 'DifOut';
actv_win   = 100;
plt_id     = 'stack_S2to3_evnt_c5';
save_fig   = 1;
fig_vis    = 'off';
fig_ftype  = 'png';
% SBJ07b_HFA_plot_stack_cond_saved_single(SBJ, elec_lab, conditions, an_id, actv_win,...
%     plt_id, save_fig, fig_vis, fig_ftype);

an_id      = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';
stat_id    = 'RL_DifOut_F01t06_WL005_WS005';%'RL_DifOut_F0t1_WL01_WS25';%
plt_id     = 'ts_F0to6_evnts_sigline';%'ts_F0to1_evnts_sigline';
for s = 1:numel(SBJs)
    SBJ08b_HFA_plot_crRT_mGLM(SBJs{s}, proc_id, an_id, stat_id, plt_id, save_fig, fig_vis, fig_ftype);
    close all;
end

%% Plot Group Stat ROI Recon
proc_id   = 'main_ft';
an_id     = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';
stat_id   = 'RL_DifOut_F01t06_WL005_WS005';%'RL_DifOut_F0t1_WL01_WS25';%
hemi      = 'l';
roi_id    = 'mgROI';
atlas_id  = 'Dx';
reg_type  = 'v';
show_lab  = 0;
save_fig  = 1;
fig_ftype = 'png';

roi_opts  = {{'l','deep',1},{'l','lat',1},{'l','MPFC',1},{'b','OFC',0}};
for roi_ix = 1:numel(roi_opts)
%     fn_view_recon_atlas_grp_stat_ROI(SBJs, proc_id, stat_id, an_id, ...
%         reg_type, show_labels, roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
%         roi_opts{roi_ix}{3},'save_fig', save_fig, 'fig_ftype', fig_ftype);
    fn_view_recon_atlas_grp_stat_venn_ROI_GLM(SBJs, proc_id, stat_id, an_id,...
        reg_type, show_lab, roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
        roi_opts{roi_ix}{3},'save_fig', save_fig, 'fig_ftype', fig_ftype);
end

%% OLDER ANOVA:
stat_id = 'corrRT_DifOutTimDO_WL200_WS50';
SBJ08a_corrRT_regressRT_ANOVA_terms(SBJ,an_id,stat_id);

%
plt_id = 'ts_S0to3_evnts_sigline';
fig_vis = 'off';
SBJ08b_HFA_plot_corrRT_ANOVA(SBJ, proc_id, an_id, stat_id, plt_id, save_fig, fig_vis, fig_ftype);

%
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


%% Plot ANOVA time series
% conditions = 'DifOut';
proc_id   = 'main_ft';
an_id     = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';
stat_id   = 'RL_DifOut_F01t06_WL005_WS005';%'RL_DifOut_F0t1_WL01_WS25';%
plt_id    = 'ts_F0to1_evnts_sigline';
atlas_id  = 'Dx';
roi_id    = 'mgROI';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';
SBJ10c_grp_GLM_plot_ts_gROIcomb(SBJs,proc_id,stat_id,an_id,...
                                    atlas_id,roi_id,plt_id,save_fig,fig_vis,fig_ftype);

%%
an_id      = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';
stat_id   = 'RL_DifOut_F01t06_WL005_WS005';%'RL_DifOut_F0t1_WL01_WS25';%
atlas_id  = 'Dx';
roi_id    = 'mgROI';
proc_id   = 'main_ft';
gm_thresh = 0;
plot_out  = 0;
plot_scat = 1;
plt_id    = '';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';
SBJ08c_HFA_grp_errbar_ROI_mGLM(SBJs,stat_id,proc_id,an_id,atlas_id,roi_id,...
                                          plot_scat,save_fig,fig_vis,fig_ftype);
              
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