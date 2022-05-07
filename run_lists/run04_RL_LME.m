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

%% Linear mixed-effects model per region
proc_id   = 'main_ft';
an_id     = 'HGm_F25t121_zbtS_sm0_l1_wn50';%'HGm_F25t121_zbtS_sm0_l1_wn100';%'HGh_F25t121_zbtS_sm0_l1';%
model_ids = {'EpnRPE_DifFB'};%, 'ERPEs_DifFB'};%{'EsRPE_DifFB','EuRPE_DifFB'};%
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
        
rcn.reg_type = 'v';
rcn.show_lab = 0;
rcn.hemi     = 'l';
rcn.mirror   = 1;

save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for m_ix = 1:numel(model_ids)
    model_id = model_ids{m_ix};
    
    for st_ix = 1:numel(stat_ids)
        stat_id = stat_ids{st_ix};
        lme_formula = lme_formulas{m_ix};
        
        % run LME
        %SBJ08g_HFA_crRT_mLME(SBJs, proc_id, an_id, model_id, stat_id, atlas_id, roi_id, lme_formula)

        % Plot fixed effects
        %SBJ08h_HFA_plot_grp_mLME(proc_id, an_id, model_id, stat_id)
        
        if strcmp(model_id, 'EpnRPE_DifFB')
            % SBJ08g_HFA_add_channel_categories(an_id, model_id, stat_id, chpval_type)
            
            % Plot ROI spider plot comparison
%               SBJ08i_HFA_plot_grp_mLME_ROI_spider(SBJ_id,proc_id,an_id,model_ids{m_ix},stat_ids{st_ix},cat_id,...
%                                                 roi_id,save_fig,'atlas_id',atlas_id,'fig_vis',fig_vis,'fig_ftype',fig_ftype);

            % Plot venn recons for INS and MPFC
            rcn.plot_roi = 'INS';
            SBJ08h_HFA_plot_grp_mLME_cat_recon(SBJ_id, proc_id, an_id, model_ids{m_ix}, stat_ids{st_ix}, cat_id,...
            	atlas_id, roi_id, rcn);
            rcn.plot_roi = 'MPFC';
            SBJ08h_HFA_plot_grp_mLME_cat_recon(SBJ_id, proc_id, an_id, model_ids{m_ix}, stat_ids{st_ix}, cat_id,...
                atlas_id, roi_id, rcn);
        end
        % Plot single channel time courses and trace plot
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
    'ERPEs_DifFB';
    'EpnRPE_DifFB';   
    };

model_names = {
    'sRPE';
    'uRPE';
    'sRPE + uRPE';
    'pRPE + nRPE'
    };

for s = 1:numel(stat_ids)
    SBJ08j_HFA_plot_grp_mLME_model_comparison(model_ids, model_names, an_id, stat_ids{s})
end

