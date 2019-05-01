%% Run scripts for SfN 2019 abstract
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath(genpath([root_dir 'PRJ_Error/scripts/']));
addpath(genpath([root_dir 'PRJ_Error/scripts/utils/']));
addpath(ft_dir);
ft_defaults

%%
SBJs = {'CP24','IR57','IR66','IR68'};%
proc_id    = 'main_ft';
an_id      = 'HGm_S_zbtS_trl2to3001_sm0_wn100_stat3';
actv_win   = 100;
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

%%
conditions = 'DifOut';
an_id = 'ERP_F_trl2to1_stat1';
for s = 1:numel(SBJs)
    SBJ06a_ERP_stats(SBJs{s},conditions,proc_id,an_id);
end

%%
conditions = 'DifOut';
plt_id     = 'stack_S2to3_evnt_c5';

SBJ07b_HFA_plot_stack_cond_saved(SBJ, conditions, an_id, actv_win, plt_id, save_fig, fig_vis, fig_ftype);

%%
for s = 3:numel(SBJs)
    SBJ07ab_HFA_actv(SBJs{s},an_id,actv_win);
end

%%
stat_id = 'corrRT_DifOutTimDO_WL200_WS50';
SBJ10a_corrRT_regressRT_ANOVA_terms(SBJ,an_id,stat_id);

%%
plt_id = 'ts_S0to3_evnts_sigline';
fig_vis = 'off';
SBJ08b_HFA_plot_corrRT_ANOVA(SBJ, proc_id, an_id, stat_id, plt_id, save_fig, fig_vis, fig_ftype);

%%
roi_id = 'gROI';
atlas_id = 'Dx';
gm_thresh = 0;
plot_out = 0;
plot_scat = 1;
plt_id      = '';
SBJ08c_HFA_GRP_summary_errbar_perc_GMlim_actv_RT_ANOVA_ROI(SBJs,stat_id,proc_id,an_id,actv_win,roi_id,...
                                                            atlas_id,gm_thresh,plt_id,plot_out,plot_scat,save_fig,fig_vis,fig_ftype);
                                                        