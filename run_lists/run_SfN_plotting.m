%% Run scripts for SfN 2019 abstract
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath(genpath([root_dir 'PRJ_Error/scripts/']));
addpath(genpath([root_dir 'PRJ_Error/scripts/utils/']));
addpath(ft_dir);
ft_defaults

%%
% behav_SBJ = {'CP22','CP23','CP24','CP25','IR57','IR60','IR62','IR63','IR65','IR66','IR67','IR68',...
%                 'IR69','IR71','IRR72','IR74','IR75','IR76','IR77','IR78','IR79','IR82','IR84'};
% eeg_behav_SBJ = {'EP06','EP07','EP08','EP10','EP11','EP14','EP15','EP16','EP17','EP18','EP19',...
%            'EEG01','EEG02','EEG03','EEG04','EEG06','EEG07','EEG08','EEG09','EEG10','EEG12'};
SBJs = {'CP24','IR57','IR68'};%

%% Print Group Accuracy
cond_lab = {'easy','hard'};
acc_cond = zeros([numel(sbj_list) 2]);
for s = 1:numel(SBJs)
    load([root_dir 'PRJ_Error/data/' SBJs{s} '/03_events/' SBJs{s} '_trl_info_final.mat'],'trl_info');
    for cond_ix = 1:2
        acc_cond(s,cond_ix) = mean(trl_info.hit(strcmp(trl_info.cond,cond_lab{cond_ix})));
    end
end
fprintf('Accuracy (n=%i):\n',numel(SBJs));
for cond_ix = 1:2
    fprintf('\t%s = %f +/- %f\n',cond_lab{cond_ix},mean(acc_cond(:,cond_ix)),std(acc_cond(:,cond_ix)));
end

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

%% Run Model (SGE...)
an_id   = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';
stat_id = 'DifOutDO_F0t1_WL1_WS25';
for s = 1:numel(SBJs)
%     SBJ08a_crRT_nANOVA(SBJs{s},an_id,stat_id);
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
proc_id   = 'main_ft';
an_id     = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';
stat_id   = 'DifOutDO_F0t1_WL1_WS25';%'corrRT_DifOutTimDO_WL200_WS50';%
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
    fn_view_recon_atlas_grp_stat_venn_ROI(SBJs, proc_id, stat_id, an_id,...
        reg_type, show_lab, roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
        roi_opts{roi_ix}{3},'save_fig', save_fig, 'fig_ftype', fig_ftype);
end

%% Plot ANOVA time series
% conditions = 'DifOut';
proc_id   = 'main_ft';
an_id     = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';%'HGm_S_zbtS_trl2to3001_sm0_wn100_stat3';%
stat_id   = 'DifOutDO_F0t1_WL1_WS25';%'corrRT_DifOutTimDO_WL200_WS50';
plt_id    = 'ts_F0to1_evnts_sigline';
atlas_id  = 'Dx';
roi_id    = 'mgROI';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype  = 'svg';
SBJ10c_grp_ANOVA_plot_ts_gROIcomb(SBJs,proc_id,stat_id,an_id,...
                                    atlas_id,roi_id,plt_id,save_fig,fig_vis,fig_ftype);

%%
an_id     = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';
stat_id   = 'DifOutDO_F0t1_WL1_WS25';
atlas_id  = 'Dx';
roi_id    = 'mgROI';
gm_thresh = 0;
plot_out  = 0;
plot_scat = 1;
plt_id    = '';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'svg';
SBJ08c_HFA_grp_errbar_ROI_nANOVA(SBJs,stat_id,proc_id,an_id,atlas_id,roi_id,...
                                          plot_scat,save_fig,fig_vis,fig_ftype);
                                      