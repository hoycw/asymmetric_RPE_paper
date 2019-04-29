%% Run scripts for SfN 2019 abstract
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath(genpath([root_dir 'PRJ_Error/scripts/']));
addpath(genpath([root_dir 'PRJ_Error/scripts/utils/']));
addpath(ft_dir);
ft_defaults

%%
SBJ        = 'IR68';
conditions = 'DifOut';
an_id      = 'HGm_S_zbtS_trl2to3001_sm0_wn100_stat3';
actv_win   = 100;
plt_id     = 'stack_S2to3_evnt_c5';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

SBJ07b_HFA_plot_stack_cond_saved(SBJ, conditions, an_id, actv_win, plt_id, save_fig, fig_vis, fig_ftype);

%%
stat_id = 'corrRT_DifOutTimDOOT_WL200_WS50';
SBJ10a_corrRT_regressRT_ANOVA_terms(SBJ,an_id,stat_id);