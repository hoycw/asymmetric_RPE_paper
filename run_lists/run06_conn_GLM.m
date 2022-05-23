%% Run scripts for RL GLM analyses
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath([root_dir 'PRJ_Error/scripts/']);
addpath(genpath([root_dir 'PRJ_Error/scripts/utils/']));
addpath(ft_dir);
ft_defaults
%%
SBJ_id = 'preproc';
SBJs = fn_load_SBJ_list(SBJ_id);
SBJs = SBJs(2:end);
% %% subject-level stats: GLM
% proc_id   = 'main_ft';
% an_id     = 'HGm_F25t121_zbtS_sm0_l1_wn50';%'HGm_F25t121_zbtS_sm0_l1_wn100';%'HGh_F25t121_zbtS_sm0_l1';%
% model_ids = {'ERPEs_DifFB'};
% stat_ids  = {'mGLM_st0t6_conn_Xcorr'}; %'mGLM_st0t10_WL05_WS25'};%
% conn_ids = {'Xcorr'}; %'MI_win','MI_0lag','corr_win','corr_0lag'};
% pmask = 1;
% nbins = 8;
% for ci = 1:numel(conn_ids)
%     for s = 3:numel(SBJs)
%         for mi = 1:numel(model_ids)
%             for si = 1:numel(stat_ids)
%                 SBJ = SBJs{s};
%                 model_id = model_ids{mi};
%                 stat_id = stat_ids{si};
%                 conn_id = conn_ids{ci};
%                 SBJ11a_HFA_conn_mGLM(SBJ, proc_id, an_id, model_id, stat_id, conn_id, pmask, nbins)
%             end
%         end
%     end
% end

%% Group-levels stats: LME
proc_id   = 'main_ft';
an_id     = 'HGm_F25t121_zbtS_sm0_l1_wn50';%'HGm_F25t121_zbtS_sm0_l1_wn100';%'HGh_F25t121_zbtS_sm0_l1';%
model_ids = {'EpnRPE_DifFB'};%;s'ERPEs_DifFB'};
stat_id = 'mLME_st0t6_WL05_WS25';
%stat_ids  = {'mGLM_st0t600_conn_MI_win'}; %'mGLM_st0t10_WL05_WS25'};%
conn_ids = {'Xcorr'};%,'corr_win','MI_win'}; %'MI_win','MI_0lag','corr_win','corr_0lag','Xcorr'};
swap_Xcorr = 1; % flip Xcorr plot
group_colors = 1;
% pmask = 1;
% nbins = 8;
lme_formula = ['y~nRPE + pRPE + EV + (1 + pRPE + nRPE + EV | sub) +', ...
                '(1+pRPE + nRPE + EV | sub:chan)'];
for ci = 1:numel(conn_ids)
    for mi = 1:numel(model_ids)
        %for si = 1:numel(stat_ids)
        model_id = model_ids{mi};
       	%stat_id = stat_ids{si};
        conn_id = conn_ids{ci};
        %SBJ11b_HFA_conn_mLME(SBJs, proc_id, an_id, model_id, conn_id, lme_formula)
        %SBJ11c_HFA_conn_peak_stats(proc_id, an_id, model_id, conn_id, swap_Xcorr)
        %SBJ11c_HFA_conn_plot_grp_mLME(proc_id, an_id, model_id, conn_id, swap_Xcorr)
        SBJ11d_HFA_conn_plot_grp_mLME_chancoef(proc_id, an_id, model_id, conn_id, group_colors, swap_Xcorr)
        %SBJ11e_HFA_conn_plot_grp_mLME_confusion(proc_id, an_id, model_id, conn_id, swap_Xcorr, stat_id)
        %end
    end
end
