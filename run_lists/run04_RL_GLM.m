%% Run scripts for RL GLM analyses
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath([root_dir 'PRJ_Error/scripts/']);
addpath(genpath([root_dir 'PRJ_Error/scripts/utils/']));
addpath(ft_dir);
ft_defaults

%%
SBJ_id = 'preproc';%'preproc_nu';%
SBJs = fn_load_SBJ_list(SBJ_id);

%% Single SBJ RL Model
proc_id   = 'main_ft';
model_ids = {'EsRPE_DifFB','EuRPE_DifFB'};%{'ERPEs_DifFB'};%'ERB_DifFB'};%'ERBuRPE_DifFB'};%
%model_ids = {'EpnRPE_DifFB'};%'RL3D_DifFB'};

fig_vis   = 'on';
save_fig  = 1;
fig_ftype = 'png';

for mdl_ix = 1:numel(model_ids)
    for s = 1:numel(SBJs)
        % Run model
        BHV02a_RL_model(SBJs{s},proc_id,model_ids{mdl_ix});

        % Plot model fit to tolerance and outcomes/accuracy
%         BHV02b_RL_model_plot(SBJs{s},proc_id,model_ids{mdl_ix},...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    
    % Plot group model fits (overlapping sigmoids without tolerance)
%     BHV02c_RL_model_plot_grp(SBJ_id,proc_id,model_ids{mdl_ix},...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
    % Plot model predicitons by condition across group
    plt_id    = 'line_cond';
    BHV02d_RL_model_plot_grp_predictions(SBJ_id,proc_id,model_ids{mdl_ix},plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    %close all;
end


%% Run Model on HFA (SGE...)
proc_id   = 'main_ft';
an_id     = 'HGm_F25t121_zbtS_sm0_l1_wn50';%'HGm_F25t121_zbtS_sm0_l1_wn100';%'HGh_F25t121_zbtS_sm0_l1';%
% model_ids = {'ERPEs_DifFB'};%'ERB_DifFB'};%'ERBuRPE_DifFB'};%
model_ids = {'EpunRPE_DifFB'};%'RL3D_DifFB'};
% stat_ids  = {'mGLM_st0t6_WL05_WS25'};%'mGLM_st0t10_WL05_WS25'};%'mGLM_EHNu_st0t6_WL05_WS25'};%
stat_ids = {'mGLM_st0t6_WL05_WS25'};%'mGLM_DifOut_st0t6_WL05_WS25'};
atlas_id  = 'Dx';

roi_id    = 'MPFCINS';%'gROI';%'main3';%'mgROI';%
plot_out  = 0;
plot_scat = 1;
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

% tbin_id     = 'cnts';
% z_thresh    = 0;
% gm_thresh   = 0;
% median_yn   = 0;

hemi      = 'l';
reg_type  = 'v';
show_lab  = 0;
mirror    = 1;
skip_reg  = 'EV';
roi_opts  = {{'l','INS',1},{'l','MPFC',1},{'l','lat',1},{'b','OFC',0}};%};%

for m_ix = 1:numel(model_ids)
    for st_ix = 1:numel(stat_ids)
        for s = 1:numel(SBJs)
            % Run Mass GLM Stats
%             SBJ08a_HFA_crRT_mGLM(SBJs{s},proc_id,an_id,model_ids{m_ix},stat_ids{st_ix});
            
            % Plot Mass GLM Results
%             plt_id    = 'ts_F2t1_evnts_sigline';
%             SBJ08b_HFA_plot_crRT_mGLM(SBJs{s}, proc_id, an_id, model_ids{m_ix}, stat_ids{st_ix},...
%                 plt_id, save_fig, 'atlas_id',atlas_id,'fig_vis',fig_vis,'fig_ftype',fig_ftype);%, 'elec_lab', {'LIN3-4'});
%             close all;
        end
        
%         % Plot bar graph showing proprotion of effects by ROI
%         SBJ08c_HFA_grp_errbar_ROI_mGLM(SBJ_id,proc_id,an_id,model_ids{m_ix},stat_ids{st_ix},...
%             roi_id,plot_scat,save_fig,'atlas_id',atlas_id,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
        % Plot histograms of betas by ROI
%         SBJ08c_HFA_plot_grp_mGLM_ROI_hist(SBJ_id,proc_id,an_id,model_ids{m_ix},stat_ids{st_ix},...
%             roi_id,save_fig,'atlas_id',atlas_id,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
        % Plot latency time series per ROI
        plt_id    = 'ts_F0t6_evnts_sigline';
        SBJ08d_HFA_plot_grp_GLM_ts_ROI_butt(SBJ_id,proc_id,an_id,model_ids{m_ix},stat_ids{st_ix},...
            roi_id,plt_id,save_fig,'atlas_id',atlas_id,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
        SBJ08d_HFA_plot_grp_GLM_ts_gROIcomb(SBJ_id,proc_id,an_id,model_ids{m_ix},stat_ids{st_ix},...
            roi_id,plt_id,save_fig,'atlas_id',atlas_id,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
        % Plot onset latencies per effect and ROI
        plt_id      = 'onsets_0t6_violin_all';
%         SBJ08f_HFA_plot_grp_GLM_onsets_ROI(SBJ_id,proc_id,an_id,model_ids{m_ix},stat_ids{st_ix},...
%             roi_id,plt_id,save_fig,'atlas_id',atlas_id,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
        % flip it and do within ROI regressor onsets!
%         SBJ08f_HFA_plot_grp_GLM_onsets_wiROI(SBJ_id,proc_id,an_id,model_ids{m_ix},stat_ids{st_ix},...
%             roi_id,plt_id,save_fig,'atlas_id',atlas_id,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
        % Recon plots
        % Group Stat Recon Viewing
        % fn_view_recon_atlas_grp_stat(SBJs, proc_id, stat_id, an_id, 'v', show_lab, 'r', atlas_id, roi_id, plot_out);
        % fn_view_recon_atlas_grp_stat(SBJs, proc_id, stat_id, an_id, 'v', show_lab, 'l', atlas_id, roi_id, plot_out);

        for roi_ix = 1:numel(roi_opts)
%                 fn_view_recon_atlas_grp_stat_ROI(SBJs, proc_id, stat_ids{st_ix}, an_id, ...
%                     reg_type, show_lab, roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
%                     roi_opts{roi_ix}{3},'save_fig', save_fig, 'fig_ftype', fig_ftype);
%             fn_view_recon_atlas_grp_stat_venn_ROI_GLM(SBJ_id, proc_id, an_id, model_ids{m_ix}, stat_ids{st_ix},...
%                 reg_type, show_lab, roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
%                 roi_opts{roi_ix}{3},'save_fig', save_fig, 'fig_ftype', fig_ftype);%, 'skip_reg', skip_reg);
            
%             % Skip EV
%             fn_view_recon_atlas_grp_stat_venn_ROI_GLM(SBJ_id, proc_id, an_id, model_ids{m_ix}, stat_ids{st_ix},...
%                 reg_type, show_lab, roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
%                 roi_opts{roi_ix}{3},'save_fig', save_fig, 'fig_ftype', fig_ftype, 'skip_reg', skip_reg);
        end
        
        % Older SBJ recon plots
        for s = 1:numel(SBJs)
            % %     fn_view_recon_stat(SBJs{s},proc_id,stat_id,an_id, 'pat', '', 0, 'r', 0, 'save_fig', 1);
            %     fn_view_recon_stat(SBJs{s},proc_id,stat_id,an_id, 'pat', '', 0, 'l', mirror, 0, 'save_fig', 1);
            %close all;
        end
    end
end

%% Plot Venn Diagram of GLM Regressor Effects by ROI
proc_id   = 'main_ft';
% stat_regs: {{'an_id','model_id','reg','stat_id'},...}
stat_regs = {{'HGm_F25t121_zbtS_sm0_l1_wn50','ERBuRPE_DifFB','uRPE','mGLM_st0t6_WL05_WS25'},...
             {'HGm_F25t121_zbtS_sm0_l1_wn50','ERBuRPE_DifFB','ERB','mGLM_st0t6_WL05_WS25'}};
% stat_regs = {{'HGm_F25t121_zbtS_sm0_l1_wn50','ERPEs_DifFB','sRPE','mGLM_st0t6_WL05_WS25'},...
%              {'HGm_F25t121_zbtS_sm0_l1_wn50','ERPEs_DifFB','uRPE','mGLM_st0t6_WL05_WS25'}};

hemi      = 'b';
atlas_id  = 'Dx';
roi_id    = 'MPFCINS';%'main3';%'gROI';
plot_out  = 0;
plt_id    = 'venn';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'svg';

SBJ08e_HFA_plot_grp_GLM_reg_pie_ROI(SBJ_id, proc_id, stat_regs, hemi, roi_id,...
                                       plot_out, plt_id, save_fig,'atlas_id',atlas_id,...
                                       'fig_vis',fig_vis,'fig_ftype',fig_ftype);

SBJ08e_HFA_plot_grp_GLM_reg_venn_ROI(SBJ_id, proc_id, stat_regs, hemi, roi_id,...
                                       plot_out, plt_id, save_fig,'atlas_id',atlas_id,...
                                       'fig_vis',fig_vis,'fig_ftype',fig_ftype);
                                   
%% ERPs: Run models and Plot Group GLM Proportion of Effects by ROI
proc_id   = 'main_ft';
% an_id     = 'ERP_F25t1';
model_ids = {'ERPEs_DifFB'};
stat_ids  = {'mGLM_st0t10_ds200'};
atlas_id  = 'Dx';
roi_id    = 'mgROI';%'main3';%'MPFCINS';%
gm_thresh = 0;
plot_out  = 0;
plot_scat = 1;

save_fig  = 1;
fig_vis   = 'off';
fig_ftype = 'png';

% tbin_id     = 'cnts';
% z_thresh    = 0;
% gm_thresh   = 0;
% median_yn   = 0;

hemi      = 'l';
reg_type  = 'v';
show_lab  = 0;
skip_reg  = 'EV';
roi_opts  = {{'l','deep',1},{'l','lat',1},{'l','MPFC',1},{'b','OFC',0}};

for m_ix = 1:numel(model_ids)
    for st_ix = 1:numel(stat_ids)
        for s = 1:numel(SBJs)
            % Run Mass GLM Stats
%             SBJ06c_ERP_crRT_mGLM(SBJs{s},proc_id,an_id,model_ids{m_ix},stat_ids{st_ix});
            
            % Plot Mass GLM Results
%             plt_id    = 'ts_F2t1_evnts_sigline';
%             SBJ06d_ERP_plot_crRT_mGLM(SBJs{s}, proc_id, an_id, model_ids{m_ix}, stat_ids{st_ix},...
%                 plt_id, save_fig, 'atlas_id',atlas_id,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%             close all;
        end
%         
%         % Plot bar graph showing proprotion of effects by ROI
%         SBJ06e_ERP_grp_errbar_ROI_mGLM(SBJ_id,proc_id,an_id,model_ids{m_ix},stat_ids{st_ix},...
%             roi_id,plot_scat,save_fig,'atlas_id',atlas_id,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%         
%         % Plot latency time series per ROI
%         plt_id    = 'ts_F0t1_evnts_sigline';
%         SBJ06f_ERP_plot_grp_GLM_ts_gROIcomb(SBJ_id,proc_id,an_id,model_ids{m_ix},stat_ids{st_ix},...
%             roi_id,plt_id,save_fig,'atlas_id',atlas_id,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%         
%         % Plot onset latencies per effect and ROI
%         plt_id      = 'onsets_0t1_violin_all';
%         SBJ06g_ERP_plot_grp_GLM_onsets_ROI(SBJ_id,proc_id,an_id,model_ids{m_ix},stat_ids{st_ix},...
%             roi_id,plt_id,save_fig,'atlas_id',atlas_id,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%         
        % flip it and do within ROI regressor onsets!
        
        % Stat Recons
        for roi_ix = 1:numel(roi_opts)
            % Plot significant electrodes
%             fn_view_recon_atlas_grp_stat_ROI(SBJs, proc_id, stat_id, an_id, ...
%                 reg_type, show_labels, roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
%                 roi_opts{roi_ix}{3},'save_fig', save_fig, 'fig_ftype', fig_ftype);
            
%             % Plot venn recons
%             fn_view_recon_atlas_grp_stat_venn_ROI_GLM(SBJ_id, proc_id, an_id, model_ids{m_ix}, stat_ids{st_ix},...
%                 reg_type, show_lab, roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
%                 roi_opts{roi_ix}{3},'save_fig', save_fig, 'fig_ftype', fig_ftype);%, 'skip_reg', skip_reg);
%             
%             % Skip EV
%             fn_view_recon_atlas_grp_stat_venn_ROI_GLM(SBJ_id, proc_id, an_id, model_ids{m_ix}, stat_ids{st_ix},...
%                 reg_type, show_lab, roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
%                 roi_opts{roi_ix}{3},'save_fig', save_fig, 'fig_ftype', fig_ftype, 'skip_reg', skip_reg);
        end
    end
end

%% Linear mixed-effects model per region
proc_id   = 'main_ft';
an_id     = 'HGm_F25t121_zbtS_sm0_l1_wn50';%'HGm_F25t121_zbtS_sm0_l1_wn100';%'HGh_F25t121_zbtS_sm0_l1';%
model_ids = {'EpnRPE_DifFB'};%,'EsRPE_DifFB','EuRPE_DifFB'};%, 'ERPEs_DifFB'};%{'EsRPE_DifFB','EuRPE_DifFB'};%
stat_ids  = {'mLME_st0t6_WL05_WS25'};%'mGLM_st0t10_WL05_WS25'};%
cat_id    = 'puns';
atlas_id  = 'Dx';
roi_id    = 'MPFCINS';%'gROI';%'mgROI';%'MPFCINS';%

lme_formulas = {['y~pRPE + nRPE + EV + (1 + pRPE + nRPE + EV | sub) +',...
                '(1+pRPE + nRPE + EV | sub:chan)']};
chpval_type = 'p'; % or q for qvals
% lme_formulas = {['y~pRPE + nRPE + EV + (1 + pRPE + nRPE + EV | sub) +',...
%     '(1+pRPE + nRPE + EV | sub:chan)'];
%     ['y~uRPE + sRPE + EV + (1 + uRPE + sRPE + EV | sub) +', ...
%     '(1+uRPE + sRPE + EV | sub:chan)']
%     };
% 
% lme_formulas = {['y~sRPE + EV + (1 + sRPE + EV | sub) +',...
%     '(1+sRPE + EV | sub:chan)'];
%     ['y~uRPE  + EV + (1 + uRPE + EV | sub) +', ...
%     '(1+uRPE + EV | sub:chan)']
%     };

%lme_formula = 'y~nRPE + pRPE + EV + (1 + pRPE + nRPE + EV | chan)';
        
reg_type  = 'v';
show_lab  = 0;
roi_opts  = {{'l','deep',1},{'l','MPFC',1}};%,{'l','lat',1},{'b','OFC',0}
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for m_ix = 1:numel(model_ids)
    model_id = model_ids{m_ix};
    
    for st_ix = 1:numel(stat_ids)
        stat_id = stat_ids{st_ix};
        %lme_formula = lme_formulas{m_ix};
        
        % run LME
        %SBJ08g_HFA_crRT_mLME(SBJs, proc_id, an_id, model_id, stat_id, atlas_id, roi_id, lme_formula)

        % Plot fixed effects
        %SBJ08h_HFA_plot_grp_mLME(proc_id, an_id, model_id, stat_id)
        
        if strcmp(model_id, 'EpnRPE_DifFB')
            %SBJ08g_HFA_add_channel_categories(an_id, model_id, stat_id, chpval_type)
            SBJ08k_HFA_channel_categories_stats(an_id, model_id, stat_id)

            % Plot ROI spider plot comparison
%               SBJ08i_HFA_plot_grp_mLME_ROI_spider(SBJ_id,proc_id,an_id,model_ids{m_ix},stat_ids{st_ix},cat_id,...
%                                                 roi_id,save_fig,'atlas_id',atlas_id,'fig_vis',fig_vis,'fig_ftype',fig_ftype);

            % Plot venn recons
    %         for roi_ix = 1:numel(roi_opts)
    %             fn_view_recon_atlas_grp_stat_ROI_LME_cat(SBJ_id, proc_id, an_id, model_ids{m_ix}, stat_ids{st_ix}, cat_id,...
    %                 reg_type, show_lab, roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
    %                 roi_opts{roi_ix}{3});
    %         end
        end
        % Plot single channel time courses, trace plot and star plots
        
        %SBJ08h_HFA_plot_grp_mLME_chancoef(proc_id, an_id, model_id, stat_id)
        % Plot example LME results in single electrode
%         SBJ08h_HFA_plot_crRT_LME('IR87', proc_id, an_id, model_ids{m_ix}, stat_ids{st_ix},...
%             'ts_F2t1_evnts_sigline', save_fig, 'atlas_id',atlas_id,'fig_vis',fig_vis,...
%             'fig_ftype',fig_ftype, 'elec_lab', {'LIN3-4'});
    end
end

%% Plot LME comparison
proc_id   = 'main_ft';
an_id     = 'HGm_F25t121_zbtS_sm0_l1_wn50';%'HGm_F25t121_zbtS_sm0_l1_wn100';%'HGh_F25t121_zbtS_sm0_l1';%
stat_ids  = {'mLME_st0t6_WL05_WS25'};%'mGLM_st0t10_WL05_WS25'};%

model_ids = {
    'EsRPE_DifFB'
    'EuRPE_DifFB';
    %'ERPEs_DifFB';
    'EpnRPE_DifFB';   
    };

model_names = {
    'Signed RPE';
    'Unsigned RPE';
    %'sRPE + uRPE';
    'Distributional RPE'
    };
colors = {[209 151 105]./255, [0 158 115]./255, [126   47  142]./255};
      
for s = 1:numel(stat_ids)
    stat_id = stat_ids{s};
    SBJ08j_HFA_plot_grp_mLME_model_comparison(model_ids, model_names, an_id, stat_id, colors)
end


