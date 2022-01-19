function SBJ08a_HFA_crRT_mLME(SBJs, proc_id, an_id, model_id, stat_id, atlas_id, roi_id)
%% function SBJ08a_crRT_mGLM(SBJ,an_id,stat_id)
%   Run GLM with sliding windows for given time-frequency analysis
%   Correlation with RT and regression of RT as confound supported
%   Model based on log fit to behavior:
%       reward prediction, signed PE, unsigned PE
% INPUTS:
%   SBJ [str] - subject ID
%   proc_id [str] - ID of processing pipeline
%   an_id [str] - ID of HFA analysis to run stats
%   model_id [str] - ID of the model used in GLM
%   stat_id [str] - ID of the statistical parameters to extract from HFA
% OUTPUTS:
%   beta [struct] - pseudo-FT structure with main GLM outputs
%   rt [FT struct] - output of correlation with RT if st.rt_corr==1
%   st [struct] - stat params loaded via stat_id

%% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

% rng('shuffle'); % seed randi with time

%% Load Data
hfa_data = [];
for s = 1:numel(SBJs)
    SBJ = SBJs{s};
    fprintf('loading subject %s\n',SBJ)
    eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
    % eval(['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m']);
    eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
    eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);
    if ~strcmp(st.an_style,'mLME'); error('This script is only for mass mLME!'); end
    
    % Load behavioral data, model, and HFA
    load([SBJ_vars.dirs.events,SBJ,'_bhv_',proc_id,'_final.mat']);
    load([SBJ_vars.dirs.models SBJ '_model_' mdl.model_id '.mat']);
    load([SBJ_vars.dirs.proc,SBJ,'_ROI_',proc_id,'_',an_id,'.mat']);
    
    % Check if more than one frequency, error for now
    if numel(hfa.freq)>1
        error('HFA has more than one frequency, can''t run on that for now!');
    end
    %% Select channels in ROIs
    [roi_list,~, roi_field] = fn_roi_label_styles(roi_id);
    load([SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_pat_' atlas_id '_final.mat']);
    r_ix = [];
    for r = 1:numel(roi_list)
        cr_ix = find(strcmp(elec.(roi_field), roi_list{r}));
        r_ix = [r_ix; cr_ix];
    end
    %% Select conditions of interest
    [reg_lab, ~, ~, ~, ~] = fn_regressor_label_styles(mdl.model_lab);
    [cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(st.stat_cond);
    [mdl_cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(mdl.model_cond);
    if ~all(contains(cond_lab,mdl_cond_lab))
        error([model_id ' does not cover all conditions in ' stat_id '!']);
    end
    full_cond_idx = fn_condition_index(cond_lab, bhv);
    bhv = fn_select_bhv(bhv, full_cond_idx);
    model = model(full_cond_idx~=0,:);
    
    % Select data in stat window
    cfg_trim = [];
    cfg_trim.latency = st.stat_lim;
    cfg_trim.trials  = full_cond_idx~=0;
    cfg_trim.channel =  elec.label(r_ix);
    hfa = ft_selectdata(cfg_trim,hfa);
    hfa_time = hfa.time;

    roi_labs = elec.(roi_field)(r_ix);
    chan_labs = elec.label(r_ix);
    %% Compute max z score per condition (for later exclusions)
    % cond_idx = fn_condition_index(cond_lab, bhv);
    %
    % max_z = zeros(size(hfa.label));
    % cfg_avg = [];
    % cfg_avg.avgoverrpt = 'yes';
    % for cond_ix = 1:numel(cond_lab)
    %     % Average within condition
    %     cfg_avg.trials = find(cond_idx==cond_ix);
    %     tmp_avg = ft_selectdata(cfg_avg,hfa);
    %     % Binarize channels by average HFA vs. threshold
    %     max_z = max(max_z,max(squeeze(tmp_avg.powspctrm),[],2));
    % end
    
    %% Run correlations with RT
    if st.rt_corr
        fprintf('================== Running Correlation with RT =======================\n');
        cfg_rt.design           = zscore(bhv.rt);
        cfg_rt.ivar             = 1;
        rt = ft_freqstatistics(cfg_rt, hfa);
    end  
    %% Regress off Reaction Time
    if st.regress_rt
        fprintf('================== Regressing RT =======================\n');
        cfg_conf = [];
        cfg_conf.output = 'residual';
        cfg_conf.confound = zscore(bhv.rt);
        hfa = ft_regressconfound(cfg_conf, hfa);
    end
    %% Average HFA in Sliding Windows
    fprintf('================== Averaging HFA within Windows =======================\n');
    win_lim = fn_get_win_lim_from_center(hfa.time,st.win_center,st.win_len);
    
    % Average in windows
    hfa_win = zeros([size(hfa.powspctrm,1) size(hfa.powspctrm,2) size(win_lim,1)]);%size(hfa.powspctrm,3)
    for w_ix = 1:size(win_lim,1)
        hfa_win(:,:,w_ix) = squeeze(nanmean(hfa.powspctrm(:,:,1,win_lim(w_ix,1):win_lim(w_ix,2)),4));
    end   
    %% Store data in table per roi
    ntimes = size(hfa_win,3);
    
    for r = 1:numel(roi_list)
        croi = roi_list{r};
        if s == 1
            hfa_data{r} = cell(numel(SBJs), ntimes);
        end
        cur_chan = chan_labs(strcmp(roi_labs,croi));
        nch = numel(cur_chan);
        for ct = 1:ntimes
            hfa_data{r}{s,ct} = cell(size(hfa_win,1)*nch,6);
            dcount = 0;
            for ch = 1:nch
                lab = cur_chan{ch};
                cdata = squeeze(hfa_win(:,ch,ct));
                for cdp = 1:numel(cdata)
                    dcount = dcount + 1;
                    hfa_data{r}{s,ct}{dcount,1} = cdata(cdp);
                    hfa_data{r}{s,ct}{dcount,2} = SBJ;
                    hfa_data{r}{s,ct}{dcount,3} = lab;
                    hfa_data{r}{s,ct}{dcount,4} = model(cdp,1);
                    hfa_data{r}{s,ct}{dcount,5} = model(cdp,2);
                    hfa_data{r}{s,ct}{dcount,6} = model(cdp,3);
                end
            end
        end
    end
    clear dirs_fields field_ix SBJ_id SBJ_vars hfa hfa_win
end

%% put together in a table
hfa_tables = [];
for r = 1:numel(hfa_data)
    hfa_tables{r} = cell(size(hfa_data{r},2),1);
    for ct = 1:size(hfa_data{r},2)
        ctable = cell2table(vertcat(hfa_data{r}{:,ct}));
        ctable.Properties.VariableNames = {'y','sub','chan','EV','sRPE','uRPE'};
        hfa_tables{r}{ct} = ctable;
    end
end
clear hfa_data

%% Run LMEs
LMEs = {};
coefs = {};
lower = {};
upper = {};
pvals = {};
lme_formula = ['y~sRPE + uRPE + EV + (1 + sRPE + uRPE + EV | sub) +', ...
                '(1 + sRPE + uRPE + EV | sub:chan)'];
for r = 1:numel(hfa_tables)
    nt = size(hfa_tables{r},1);
    LMEs{r} = cell(nt,1);
    pvals{r} = NaN(size(model,2) + 1,nt);
    coefs{r} = NaN(size(model,2) + 1,nt);
    lower{r} = NaN(size(model,2) + 1,nt);
    upper{r} = NaN(size(model,2) + 1,nt);
    for t = 1:nt
        fprintf('fitting label %d and time point %d\n',r,t)
        LMEs{r}{t} = fitlme(hfa_tables{r}{t}, lme_formula);
        pvals{r}(:,t) = LMEs{r}{t}.Coefficients.pValue;
        coefs{r}(:,t) = LMEs{r}{t}.Coefficients.Estimate;
        lower{r}(:,t) = LMEs{r}{t}.Coefficients.Lower;
        upper{r}(:,t) = LMEs{r}{t}.Coefficients.Upper;
    end
end
%% FDR correction
qvals = cell(size(pvals));
for r = 1:numel(qvals)
    qvals{r} = NaN(size(pvals{r}));
    qvals{r}(1,:) = pvals{r}(1,:); % don't include intercept in correction
    for cf = 1:size(pvals{r},2)
        [~, ~, ~, qvals{r}(2:end,:)] = fdr_bh(pvals{r}(2:end,:));
    end
end
%% Store in a structure
beta =  [];
%conn_stats.LMEs = LMEs;
beta.feature   = reg_lab;
beta.time      = st.win_center;
beta.coefs = coefs;
beta.lower = lower;
beta.upper = upper;
beta.pvals = pvals;
beta.qvals = qvals;
beta.label = roi_list;
beta.win_lim   = win_lim;
beta.win_lim_s = hfa_time(win_lim);

%% Save
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
out_fname = [stats_dir model_id '_' stat_id '_' an_id '_hfa.mat'];
save(out_fname,'-v7.3','beta')

%% Run ANOVA
% fprintf('================== Running Mass LME =======================\n');
% % Create structure for beta in fieldtrip style
% if ~any(strcmp(reg_lab,'offset'))
%     reg_lab = ['offset' reg_lab];
%     model   = [ones(size(model,1), 1) model];
% end
% beta.model     = model;
% beta.feature   = reg_lab;
% beta.time      = st.win_center;
% beta.win_lim   = win_lim;
% beta.label     = hfa.label;
% beta.dimord    = 'rpt_chan_time';
% beta.max_hfa_z = max_z;
% % w2.trial     = zeros([numel(w2.cond) length(hfa.label) length(w2.time)]);
% beta.boot      = zeros([numel(beta.feature) length(hfa.label) length(beta.time) st.n_boots]);
% % w2.pval      = w2.trial;
% beta.win_lim_s = hfa.time(win_lim);
% 
% % Compute GLM for real model
% beta.trial = fn_mass_GLM(beta.model,hfa_win,0);
% 
% % Compute GLM for permuted data
% % b = '';
% fprintf('boot #: ');
% for boot_ix = 1:st.n_boots
% %     m = sprintf(' permutation %d/%d', boot_ix, n_boots);
% %     fprintf([b m]); b = repmat('\b',[1 length(m)]);
%     fprintf('%i..',boot_ix);
%     rand_model = model(randperm(size(model,1)),:);
%     beta.boot(:,:,:,boot_ix) = fn_mass_GLM(rand_model,hfa_win,0);
%     if mod(boot_ix,20)==0
%         fprintf('\n');
%     end
% end
% 
% % Compute statistics
% % Original one-sided logic (e.g., ANOVA variance explained): % sum(boots>real)/n_boots
% %   one-sided code: beta.pval = sum(bsxfun(@ge,beta.boot,beta.trial),4)./st.n_boots;
% % Here, using two-sided logic: sum(abs(boots-boot_median)>abs(real-boot_median))/n_boots
% %   Note the bootstrapped distribution and real beta are median-centered to
% %   correct for bias before mirroring the distribution to allow one-sided
% %   FDR correction and testing
% beta.pval = sum(bsxfun(@ge,abs(beta.boot-median(beta.boot,4)),abs(beta.trial-median(beta.boot,4))),4)./st.n_boots;
% 
% % Multiple Comparisons Correction within Channel
% beta.qval = nan(size(beta.pval));
% for ch_ix = 1:numel(beta.label)
%     [~, ~, ~, beta.qval(:,ch_ix,:)] = fdr_bh(squeeze(beta.pval(:,ch_ix,:)));
% end
% 
% % Z-score beta statistics for later group use
% beta.zscore   = norminv(1-beta.pval,0,1);
% beta.bootmean = mean(beta.boot,4);
% beta.bootstd  = std(beta.boot,[],4);
% beta.zscore(isinf(beta.zscore)) = norminv(1-1/st.n_boots/2,0,1);
% 
% %% Print results
% % Prep report
% sig_report_fname = [SBJ_vars.dirs.stats SBJ '_mGLM_ROI_' model_id '_' stat_id '_' an_id '_sig_report.txt'];
% if exist(sig_report_fname)
%     system(['mv ' sig_report_fname ' ' sig_report_fname(1:end-4) '_bck.txt']);
% end
% sig_report = fopen(sig_report_fname,'a');
% if st.rt_corr; reg_lab = [reg_lab {'RT'}]; end
% result_str = ['%-10s' repmat('%-10i',[1 numel(reg_lab)]) '\n'];
% 
% % Print header
% fprintf(sig_report,'%s (n = %i)\n',SBJ,numel(beta.label));
% fprintf(sig_report,[repmat('%-10s',[1 1+numel(reg_lab)]) '\n'],'label',reg_lab{:});
% 
% % Print summary lines (absolute)
% if st.rt_corr
%     fprintf(sig_report,result_str, 'count', sum(any(beta.qval(:,:,:)<st.alpha,3),2), sum(any(rt.mask(:,1,:),3)));
%     fprintf(sig_report,strrep(result_str,'i','.3f'), 'percent',...
%         [sum(any(beta.qval(:,:,:)<st.alpha,3),2)', sum(any(rt.mask(:,1,:),3))]./numel(beta.label));
% else
%     fprintf(sig_report,result_str, 'count', sum(any(beta.qval(:,:,:)<st.alpha,3),2));
%     fprintf(sig_report,strrep(result_str,'i','.3f'), 'percent',...
%         sum(any(beta.qval(:,:,:)<st.alpha,3),2)'./numel(beta.label));
% end
% 
% % Print Channel Lines
% sig_mat = zeros([numel(beta.label) numel(reg_lab)]);
% for ch_ix = 1:numel(beta.label)
%     % Consolidate to binary sig/non-sig
%     for reg_ix = 1:numel(reg_lab)
%         if any(squeeze(beta.qval(reg_ix,ch_ix,:))<st.alpha)
%             sig_mat(ch_ix,reg_ix) = 1;
%         end
%     end
%     if st.rt_corr && any(rt.mask(ch_ix,1,:))
%         sig_mat(ch_ix,numel(reg_lab)+1) = 1;
%     end
%     
%     % Report on significant electrodes for this SBJ
%     fprintf(sig_report,result_str,beta.label{ch_ix},sig_mat(ch_ix,:));
% end
% 
% fclose(sig_report);
% 
% %% Save Results
% out_fname = [SBJ_vars.dirs.stats SBJ '_mGLM_ROI_' model_id '_' stat_id '_' an_id '.mat'];
% if st.rt_corr
%     save(out_fname,'-v7.3','beta','rt');
% else
%     save(out_fname,'-v7.3','beta');
% end

end