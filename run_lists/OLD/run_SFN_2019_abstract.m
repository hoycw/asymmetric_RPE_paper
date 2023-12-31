%% Run scripts for SfN 2019 abstract
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath(genpath([root_dir 'PRJ_Error/scripts/']));
addpath(genpath([root_dir 'PRJ_Error/scripts/utils/']));
addpath(ft_dir);
ft_defaults

%%
SBJs = {'CP24','IR57','IR66','IR68'};
proc_id    = 'main_ft';
an_id      = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';%'HGm_S_zbtS_trl2to3001_sm0_wn100_stat3';
actv_win   = 100;
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

%%
an_id = 'ERP_F_trl2to1_stat1';
for s = 2:numel(SBJs)
    conditions = 'EzOut';
    SBJ06a_ERP_stats(SBJs{s},conditions,proc_id,an_id);
    conditions = 'HdOut';
    SBJ06a_ERP_stats(SBJs{s},conditions,proc_id,an_id);
end

%%
fig_vis = 'off';
conditions = 'DifOut';
an_id = 'ERP_F_trl2to1_stat1';
for s = 1:numel(SBJs)
    plt_id     = 'ts_trl2to1000_errbr';
    SBJ06b_ERP_plot_stats(SBJs{s},conditions,proc_id,an_id,plt_id,save_fig,fig_vis)
    close all;
end

%%
for s = 1:numel(SBJs)
    SBJ07ab_HFA_actv(SBJs{s},an_id,actv_win)
end

%%
fig_vis = 'off';
conditions = 'DifOut';
an_id      = 'HGm_S_zbtS_trl2to3001_sm0_wn100_stat3';%'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';%
stat_id = 'corrRT_DifOutTimDO_WL200_WS50';
for s = 1:numel(SBJs)
%     plt_id     = 'stack_S2to3_evnt_c5';%'stack_F2to12_evnt_c5';
%     SBJ07b_HFA_plot_stack_cond_saved(SBJs{s}, conditions, an_id, actv_win, plt_id, save_fig, fig_vis, fig_ftype);
%     close all;
    
    plt_id = 'ts_S0to3_evnts_sigline';%'ts_F0to1_evnts_sigline';
    SBJ08b_HFA_plot_corrRT_ANOVA(SBJs{s}, proc_id, an_id, stat_id, plt_id, save_fig, fig_vis, fig_ftype);
    close all;
end

%%
for s = 1:numel(SBJs)
    SBJ07ab_HFA_actv(SBJs{s},an_id,actv_win);
end

%%
stat_id = 'corrRT_DifOutTimDO_WL200_WS50';
SBJ08a_corrRT_regressRT_ANOVA_terms(SBJ,an_id,stat_id);

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
an_id      = 'HGm_S_zbtS_trl2to3001_sm0_wn100_stat3';%'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';%
s = 4;
fn_view_recon_stat(SBJs{s},proc_id,stat_id,an_id, 'pat', '', 0, 'r', 0);
fn_view_recon_stat(SBJs{s},proc_id,stat_id,an_id, 'pat', '', 0, 'l', 0);

%%
SBJs = {'CP24','IR57','IR68'};%'IR66' - doesn't have mni_v_Dx yet
roi_id = 'gROI';
atlas_id = 'Dx';
plot_out = 0;
show_lab = 0;
fn_view_recon_atlas_grp_stat(SBJs, proc_id, stat_id, an_id, 'v', show_lab, 'r', atlas_id, roi_id, plot_out);
fn_view_recon_atlas_grp_stat(SBJs, proc_id, stat_id, an_id, 'v', show_lab, 'l', atlas_id, roi_id, plot_out);