%% Run scripts for SfN 2019 abstract
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath(genpath([root_dir 'PRJ_Error/scripts/']));
addpath(genpath([root_dir 'PRJ_Error/scripts/utils/']));
addpath(ft_dir);
ft_defaults

%%
SBJ_id = 'preproc';
SBJs = fn_load_SBJ_list(SBJ_id);

%% Single SBJ RL Model
proc_id   = 'main_ft';
model_ids = {'ERPEs_DifFB'};

fig_vis   = 'on';
save_fig  = 1;
fig_ftype = 'png';

for mdl_ix = 1:numel(model_ids)
    for s = 1:numel(SBJs)
        % Run model
%         BHV02a_RL_model(SBJs{s},proc_id,model_ids{mdl_ix});
% 
%         % Fig. 1D: Plot model fit to tolerance and outcomes/accuracy
%         BHV02b_RL_model_plot(SBJs{s},proc_id,model_ids{mdl_ix},...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
%     
%     % Fig. 1D inset: Plot group model fits (overlapping sigmoids without tolerance)
%     BHV02c_RL_model_plot_grp(SBJ_id,proc_id,model_ids{mdl_ix},...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%     
%     % Sup. Fig. 1: Plot model predicitons by condition across group
%     plt_id    = 'line_cond';
%     BHV02d_RL_model_plot_grp_predictions(SBJ_id,proc_id,model_ids{mdl_ix},plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%     %close all;
end


%% Run Model (SGE...)
proc_id   = 'main_ft';
an_id     = 'HGm_F25t121_zbtS_sm0_l1_wn100';%'HGh_F25t121_zbtS_sm0_l1';%
model_ids = {'ERPEs_DifFB'};
stat_ids  = {'mGLM_st05t6_WL05_WS25','mGLM_st05t10_WL1_WS25'};
atlas_id  = 'Dx';
% stat_ids  = {'lme_st05t5'};%'ERPEsL_DifFB_lme_st05t5'};
% stat_id = 'RL_DifOut_F0t1_WL01_WS25';%'ERPEs_DO_glm_F0t5_WL01_WS25';

plt_id    = 'ts_F0t1_evnts_sigline';
save_fig  = 1;
fig_vis   = 'off';
fig_ftype = 'png';

for m_ix = 1:numel(model_ids)
    for st_ix = 1:numel(stat_ids)
        for s = 1:numel(SBJs)
            % Run Mass GLM Stats
%             SBJ08a_crRT_mGLM(SBJs{s},proc_id,an_id,model_ids{m_ix},stat_ids{st_ix});
            
            % Plot Mass GLM Results
            SBJ08b_HFA_plot_crRT_mGLM(SBJs{s}, proc_id, an_id, model_ids{m_ix}, stat_ids{st_ix},...
                plt_id, save_fig, 'atlas_id',atlas_id,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
            close all;
        end
    end
end

%% Plot Group GLM Proportion of Effects by ROI
proc_id   = 'main_ft';
an_id     = 'HGm_F25t121_zbtS_sm0_l1_wn100';%'HGh_F25t121_zbtS_sm0_l1';%
model_ids = {'ERPEs_DifFB'};
stat_ids  = {'mGLM_st05t6_WL05_WS25'};%'mGLM_st05t10_WL1_WS25'};
atlas_id  = 'Dx';
roi_id    = 'mgROI';
gm_thresh = 0;
plot_out  = 0;
plot_scat = 1;
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for m_ix = 1:numel(model_ids)
    for st_ix = 1:numel(stat_ids)
        % Plot bar graph showing proprotion of effects by ROI
        SBJ08c_HFA_grp_errbar_ROI_mGLM(SBJ_id,proc_id,an_id,model_ids{m_ix},stat_ids{st_ix},...
            roi_id,plot_scat,save_fig,'atlas_id',atlas_id,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
        % Plot latency time series per ROI
        plt_id    = 'ts_F0t1_evnts_sigline';
        SBJ08d_HFA_plot_grp_GLM_ts_gROIcomb(SBJ_id,proc_id,an_id,model_ids{m_ix},stat_ids{st_ix},...
            roi_id,plt_id,save_fig,'atlas_id',atlas_id,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
end

%% Plot Group Stat ROI Recon
proc_id   = 'main_ft';
an_id     = 'HGm_F25t121_zbtS_sm0_l1_wn100';%'HGh_F25t121_zbtS_sm0_l1';%
model_ids = {'ERPEs_DifFB'};
stat_ids  = {'mGLM_st05t6_WL05_WS25'};%'mGLM_st05t10_WL1_WS25'};
hemi      = 'l';
roi_id    = 'mgROI';
atlas_id  = 'Dx';
reg_type  = 'v';
show_lab  = 0;
save_fig  = 1;
fig_ftype = 'png';

roi_opts  = {{'l','deep',1},{'l','lat',1},{'l','MPFC',1},{'b','OFC',0}};
for m_ix = 1:numel(model_ids)
    for st_ix = 1:numel(stat_ids)
        for roi_ix = 1:numel(roi_opts)
            %     fn_view_recon_atlas_grp_stat_ROI(SBJs, proc_id, stat_id, an_id, ...
            %         reg_type, show_labels, roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
            %         roi_opts{roi_ix}{3},'save_fig', save_fig, 'fig_ftype', fig_ftype);
            fn_view_recon_atlas_grp_stat_venn_ROI_GLM(SBJ_id, proc_id, an_id, model_ids{m_ix}, stat_ids{st_ix},...
                reg_type, show_lab, roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
                roi_opts{roi_ix}{3},'save_fig', save_fig, 'fig_ftype', fig_ftype);
        end
    end
end

%% Single SBJ and Group Stat Recon viewing
an_id = 'HGm_F_zbtS_trl2to1201_sm0_wn100_stat1';%'HGm_S_zbtS_trl2to3001_sm0_wn100_stat3';%
roi_id = 'gROI';
atlas_id = 'Dx';
plot_out = 0;
show_lab = 0;
mirror = 1;

for s = 1:numel(SBJs)
% %     fn_view_recon_stat(SBJs{s},proc_id,stat_id,an_id, 'pat', '', 0, 'r', 0, 'save_fig', 1);
%     fn_view_recon_stat(SBJs{s},proc_id,stat_id,an_id, 'pat', '', 0, 'l', mirror, 0, 'save_fig', 1);
    %close all;
end

% Group Stat Recon Viewing
% fn_view_recon_atlas_grp_stat(SBJs, proc_id, stat_id, an_id, 'v', show_lab, 'r', atlas_id, roi_id, plot_out);
% fn_view_recon_atlas_grp_stat(SBJs, proc_id, stat_id, an_id, 'v', show_lab, 'l', atlas_id, roi_id, plot_out);    

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

%% OLDER ANOVA:
% stat_id = 'corrRT_DifOutTimDO_WL200_WS50';
% SBJ08a_corrRT_regressRT_ANOVA_terms(SBJ,an_id,stat_id);
% 
% %
% plt_id = 'ts_S0to3_evnts_sigline';
% fig_vis = 'off';
% SBJ08b_HFA_plot_corrRT_ANOVA(SBJ, proc_id, an_id, stat_id, plt_id, save_fig, fig_vis, fig_ftype);
% 
% %
% roi_id = 'gROI';
% atlas_id = 'Dx';
% gm_thresh = 0;
% plot_out = 0;
% plot_scat = 1;
% plt_id      = '';
% fig_ftype = 'png';
% SBJ08c_HFA_GRP_summary_errbar_perc_GMlim_actv_RT_ANOVA_ROI(SBJs,stat_id,proc_id,an_id,actv_win,roi_id,...
%                                                             atlas_id,gm_thresh,plt_id,plot_out,plot_scat,save_fig,fig_vis,fig_ftype);
% SBJ08c_HFA_GRP_summary_errbar_perc_GMlim_actv_RT_ANOVA_ROI_DO(SBJs,stat_id,proc_id,an_id,actv_win,roi_id,...
%                                                             atlas_id,gm_thresh,plt_id,plot_out,plot_scat,save_fig,fig_vis,fig_ftype);

