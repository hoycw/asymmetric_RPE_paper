function SBJ08a_crRT_mGLM(SBJ,an_id,stat_id)
%% function SBJ08a_crRT_mGLM(SBJ,an_id,stat_id)
%   Run GLM with sliding windows for given time-frequency analysis
%   Correlation with RT and regression of RT as confound supported
%   Model based on log fit to behavior:
%       reward prediction, signed PE, unsigned PE
% INPUTS:
%   SBJ [str] - subject ID
%   an_id [str] - HFA analysis to run stats
%   stat_id [str] - ID of the statistical parameters and design
% OUTPUTS:
%   w2 [struct] - pseudo-FT structure with main ANOVA output
%   rt [FT struct] - output of correlation with RT if st.rt_corr==1
%   st [struct] - stat params loaded via stat_id

%% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

rng('shuffle'); % seed randi with time

%% Load Data
eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
% eval(['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);

load(strcat(SBJ_vars.dirs.events,SBJ,'_trl_info_final.mat'));
load(strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'.mat'));
[cond_lab, ~, ~, ~] = fn_condition_label_styles(st.trial_cond{1});
cond_idx = fn_condition_index(st.trial_cond{cond_ix}, trl_info);

% Check if more than one frequency, error for now
if numel(hfa.freq)>1
    error('HFA has more than one frequency, can''t run on that for now!');
end

% Select data in stat window
cfg_trim = [];
cfg_trim.latency = st.stat_lim;
hfa = ft_selectdata(cfg_trim,hfa);

%% Exclude trials based on RT and trial type
% Select trial types of interest
if ~strcmp(st.trial_cond,'all')
    good_cond_idx = cond_idx~=0;
else
    good_cond_idx = true(size(trl_info.trl_n));
end

% Log combined bad trial types
% good_trl_idx = all([~bad_rt_idx good_cond_idx],2);
% st.bad_trials.all  = trl_info.trial_n(~good_trl_idx);
% %st.bad_trials.rt   = trl_info.trial_n(bad_rt_idx);
% st.bad_trials.cond = trl_info.trial_n(~good_cond_idx);

% Exclude bad trials
ti_fields = fieldnames(trl_info);
orig_n_trials = numel(trl_info.trl_n);
for f_ix = 1:numel(ti_fields)
    if numel(trl_info.(ti_fields{f_ix}))==orig_n_trials
        trl_info.(ti_fields{f_ix}) = trl_info.(ti_fields{f_ix})(good_cond_idx);
    end
end

%% Create RL Model: Fit Behavior, Compute Regressors
% Run Logistic Regression for Win Prediction
pWin = nan(size(trl_info.trl_n));
sPE  = nan(size(trl_info.trl_n));
uPE  = nan(size(trl_info.trl_n));

% Select Data (fit on everything except surprise since no outcome)
% s_idx = fn_condition_index({'Su'},bhvs{s});
X = trl_info.tol; %(~s_idx);
y = double(trl_info.hit); %(~s_idx));

% Logistic regression
rl_betas = glmfit(X,y,'binomial','link','logit');

z = rl_betas(1) + (trl_info.tol * rl_betas(2));
pWin = 1 ./ (1+exp(-z));
expected_score = pWin*2 - 1;
sPE = double(trl_info.score)/100 - expected_score;
uPE = abs(sPE);

% Build full model
model = nan([numel(trl_info.trl_n) numel(st.regressors)]);
for reg_ix = 1:numel(st.regressors)
    if strcmp(st.regressors{reg_ix},'offset')
        model(:,reg_ix) = ones(size(trl_info.trl_n));
    else
        eval(['model(:,reg_ix) = ' st.regressors{reg_ix} ';']);
    end
end

%% Compute max z score per condition (for later exclusions)
max_z = zeros(size(hfa.label));
cfg_avg = [];
cfg_avg.avgoverrpt = 'yes';
for cond_ix = 1:numel(cond_lab)
    % Average within condition
    cfg_avg.trials = find(cond_idx==cond_ix);
    tmp_avg = ft_selectdata(cfg_avg,hfa);
    % Binarize channels by average HFA vs. threshold
    max_z = max(max_z,max(squeeze(tmp_avg.powspctrm),[],2));
end

%% Run correlations with RT
if st.rt_corr
    fprintf('================== Running Correlation with RT =======================\n');
    cfg_rt.design           = zscore(trl_info.rt);
    cfg_rt.ivar             = 1;
    rt = ft_freqstatistics(cfg_rt, hfa);
end
% OUTPUT:
%   .rho   = correlation coefficients
%   .stat  = t values
%   .prob  = p values

%% Regress off Reaction Time
if st.regress_rt
    fprintf('================== Regressing RT =======================\n');
    cfg_conf = [];
    cfg_conf.output = 'residual';
    cfg_conf.confound = zscore(trl_info.rt);
    hfa = ft_regressconfound(cfg_conf, hfa);
end
% OUTPUT:
%   As of August 2017, can no longer get residuals, model, and beta; now choose via .output argument
%   'beta'  = weights
%       .beta = 4D double
%   'model' = confounds * weights = X * X\Y
%       output is a .model struct that contains .powspctrm, .time, .trialinfo, etc.
%   'residual' = Yclean = Y - X * X\Y (i.e., the residuals after subtracting the model off
%       output is .powspctrm as normal
% NOTE: DO NOT TRUST .stat and .prob output if .statistics = 'yes'!!! (e.g., p values of 2...)

%% Average HFA in Sliding Windows
fprintf('================== Averaging HFA within Windows =======================\n');
% Sliding window parameters
win_lim    = fn_sliding_window_lim(squeeze(hfa.powspctrm(1,1,1,:)),...
    round(st.win_len*trl_info.sample_rate),...
    round(st.win_step*trl_info.sample_rate));
win_center = round(mean(win_lim,2));

% Average in windows
hfa_win = zeros([size(hfa.powspctrm,1) size(hfa.powspctrm,2) size(win_lim,1)]);%size(hfa.powspctrm,3)
for w_ix = 1:size(win_lim,1)
    hfa_win(:,:,w_ix) = squeeze(nanmean(hfa.powspctrm(:,:,1,win_lim(w_ix,1):win_lim(w_ix,2)),4));
end

%% Run ANOVA
fprintf('================== Running ANOVA =======================\n');
% Create structure for w2 in fieldtrip style
beta.model     = model;
beta.feature   = st.regressors;
beta.time      = hfa.time(win_center);
beta.win_lim   = win_lim;
beta.label     = hfa.label;
beta.dimord    = 'rpt_chan_time';
beta.max_hfa_z = max_z;
% w2.trial     = zeros([numel(w2.cond) length(hfa.label) length(w2.time)]);
beta.boot      = zeros([numel(beta.feature) length(hfa.label) length(beta.time) st.n_boots]);
% w2.pval      = w2.trial;
beta.win_lim_s = hfa.time(win_lim);

% Compute ANOVA and Explained Variance for real model
if any(strcmp(st.regressors,'offset'))
    add_offset = false;
else
    add_offset = true;
end
beta.trial = fn_mass_GLM(beta.model,hfa_win,add_offset);

% Compute ANOVA for permuted data
rand_model = beta.model;
% b = '';
fprintf('boot #: ');
for boot_ix = 1:st.n_boots
%     m = sprintf(' permutation %d/%d', boot_ix, n_boots);
%     fprintf([b m]); b = repmat('\b',[1 length(m)]);
    fprintf('%i..',boot_ix);
    rand_model = rand_model(randperm(size(rand_model,1)),:);
    beta.boot(:,:,:,boot_ix) = fn_mass_GLM(rand_model,hfa_win,add_offset);
    if mod(boot_ix,20)==0
        fprintf('\n');
    end
end

% Compute statistics
beta.pval = sum(bsxfun(@ge,beta.boot,beta.trial),4)./st.n_boots; % sum(boots>real)/n_boots
beta.zscore   = norminv(1-beta.pval,0,1);
beta.bootmean = mean(beta.boot,4);
beta.bootstd  = std(beta.boot,[],4);
% w2 = rmfield(w2,'boot');
beta.zscore(isinf(beta.zscore)) = norminv(1-1/st.n_boots/2,0,1);

% Multiple Comparisons Correction within Channel
beta.qval = nan(size(beta.pval));
for ch_ix = 1:numel(beta.label)
    [~, ~, ~, beta.qval(:,ch_ix,:)] = fdr_bh(squeeze(beta.pval(:,ch_ix,:)));
end

%% Print results
% Prep report
sig_report_fname = [SBJ_vars.dirs.proc SBJ '_mGLM_ROI_' stat_id '_' an_id '_sig_report.txt'];
if exist(sig_report_fname)
    system(['mv ' sig_report_fname ' ' sig_report_fname(1:end-4) '_bck.txt']);
end
sig_report = fopen(sig_report_fname,'a');
if st.rt_corr; cond_lab = [st.regressors {'RT'}]; else cond_lab = st.regressors; end
result_str = ['%-10s' repmat('%-10i',[1 numel(cond_lab)]) '\n'];

% Print header
fprintf(sig_report,'%s (n = %i)\n',SBJ,numel(beta.label));
fprintf(sig_report,[repmat('%-10s',[1 1+numel(cond_lab)]) '\n'],'label',cond_lab{:});

% Print summary lines (absolute)
if st.rt_corr
    fprintf(sig_report,result_str, 'count', sum(any(beta.qval(:,:,:)<st.alpha,3),2), sum(any(rt.mask(:,1,:),3)));
    fprintf(sig_report,strrep(result_str,'i','.3f'), 'percent',...
        [sum(any(beta.qval(:,:,:)<st.alpha,3),2)', sum(any(rt.mask(:,1,:),3))]./numel(beta.label));
else
    fprintf(sig_report,result_str, 'count', sum(any(beta.qval(:,:,:)<st.alpha,3),2));
    fprintf(sig_report,strrep(result_str,'i','.3f'), 'percent',...
        sum(any(beta.qval(:,:,:)<st.alpha,3),2)'./numel(beta.label));
end

% Print Channel Lines
sig_mat = zeros([numel(beta.label) numel(cond_lab)]);
for ch_ix = 1:numel(beta.label)
    % Consolidate to binary sig/non-sig
    for grp_ix = 1:numel(st.regressors)
        if any(squeeze(beta.qval(grp_ix,ch_ix,:))<st.alpha)
            sig_mat(ch_ix,grp_ix) = 1;
        end
    end
    if st.rt_corr && any(rt.mask(ch_ix,1,:))
        sig_mat(ch_ix,grp_ix) = 1;
    end
    
    % Report on significant electrodes for this SBJ
    fprintf(sig_report,result_str,beta.label{ch_ix},sig_mat(ch_ix,:));
end

fclose(sig_report);

%% Save Results
out_fname = [SBJ_vars.dirs.proc SBJ '_mGLM_ROI_' stat_id '_' an_id '.mat'];
if st.rt_corr
    save(out_fname,'-v7.3','beta','rt','st');
else
    save(out_fname,'-v7.3','beta','st');
end

end
