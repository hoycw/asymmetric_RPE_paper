function SBJ06c_ERP_crRT_mGLM(SBJ, proc_id, an_id, model_id, stat_id)
%% function SBJ06c_ERP_crRT_mGLM(SBJ,an_id,stat_id)
%   Run GLM with sliding windows for given ERP analysis
%   NOT SUPPORTED: Correlation with RT
%       regression of RT as confound is still supported
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
eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
% eval(['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);
if ~strcmp(st.an_style,'mGLM'); error('This script is only for mass GLM!'); end

% Load behavioral data, model, and HFA
load([SBJ_vars.dirs.events,SBJ,'_bhv_',proc_id,'_final.mat']);
load([SBJ_vars.dirs.models SBJ '_model_' mdl.model_id '.mat']);
load([SBJ_vars.dirs.proc,SBJ,'_ROI_',proc_id,'_',an_id,'.mat']);

%% Select conditions of interest
[reg_lab, ~, ~, ~, ~] = fn_regressor_label_styles(mdl.model_lab);
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(st.stat_cond);
[mdl_cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(mdl.model_cond);
if ~all(contains(cond_lab,mdl_cond_lab))
    error([model_id ' does not cover all conditions in ' stat_id '!']);
end
cond_idx = fn_condition_index(cond_lab, bhv);
good_cond_ix = unique(cond_idx);
bhv = fn_select_bhv(bhv, cond_idx);

% Select data in stat window
cfg_trim = [];
cfg_trim.latency = st.stat_lim;
cfg_trim.trials  = cond_idx~=0;
erp_trl = ft_selectdata(cfg_trim,erp_trl);

% Log combined bad trial types
% good_trl_idx = all([~bad_rt_idx good_cond_idx],2);
% st.bad_trials.all  = bhv.trial_n(~good_trl_idx);
% %st.bad_trials.rt   = bhv.trial_n(bad_rt_idx);
% st.bad_trials.cond = bhv.trial_n(~good_cond_idx);

%% Compute max ERP per condition (for later exclusions)
% Average ERPs to get max values
erp     = cell(size(cond_lab));
max_abs = zeros(size(erp_trl.label));
for cond_ix = 1:numel(cond_lab)
    if any(good_cond_ix==cond_ix)
        cfgavg.trials = find(cond_idx==cond_ix);
        erp{cond_ix}  = ft_timelockanalysis([],erp_trl);
        % Binarize channels by average HFA vs. threshold
        max_abs = max(max_abs,max(abs(erp{cond_ix}.avg),[],2));
    end
end

%% Run correlations with RT
if st.rt_corr
    fprintf('================== Running Correlation with RT =======================\n');
    error('check this, needs adaptation to ft_timelockstatistics');
    cfg_rt.design           = zscore(bhv.rt);
    cfg_rt.ivar             = 1;
    rt = ft_timelockstatistics(cfg_rt, erp_trl);
end
% OUTPUT:
%   .rho   = correlation coefficients
%   .stat  = t values
%   .prob  = p values

%% Regress off Reaction Time
if st.regress_rt
    fprintf('================== Regressing RT =======================\n');
    error('make sure this works before running!');
    cfg_conf = [];
    cfg_conf.output = 'residual';
    cfg_conf.confound = zscore(bhv.rt);
    erp_trl = ft_regressconfound(cfg_conf, erp_trl);
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

%% Donwsample ERPs
if st.downsample
    fprintf('================== Downsampling ERPs =======================\n');
    cfgds = [];
    cfgds.resamplefs = st.dsamp_freq;
    cfgds.detrend = 'no';
    erp_trl = ft_resampledata(cfgds, erp_trl);
end

%% Run ANOVA
fprintf('================== Running Mass GLM =======================\n');
% Create structure for beta in fieldtrip style
if ~any(strcmp(reg_lab,'offset'))
    reg_lab = ['offset' reg_lab];
    model   = [ones(size(model,1), 1) model];
end
beta.model     = model;
beta.feature   = reg_lab;
beta.time      = erp_trl.time{1};
beta.label     = erp_trl.label;
beta.dimord    = 'rpt_chan_time';
beta.max_abs   = max_abs;
beta.boot      = zeros([numel(beta.feature) length(erp_trl.label) length(beta.time) st.n_boots]);

% Extract trial data into matrix
erp_data = zeros([numel(erp_trl.trial) numel(erp_trl.label) numel(erp_trl.time{1})]);
for trl_ix = 1:numel(erp_trl.trial)
    erp_data(trl_ix,:,:) = erp_trl.trial{trl_ix};
end

% Compute ANOVA and Explained Variance for real model
beta.trial = fn_mass_GLM(beta.model,erp_data,0);

% Compute ANOVA for permuted data
rand_model = beta.model;
% b = '';
fprintf('boot #: ');
for boot_ix = 1:st.n_boots
%     m = sprintf(' permutation %d/%d', boot_ix, n_boots);
%     fprintf([b m]); b = repmat('\b',[1 length(m)]);
    fprintf('%i..',boot_ix);
    rand_model = rand_model(randperm(size(rand_model,1)),:);
    beta.boot(:,:,:,boot_ix) = fn_mass_GLM(rand_model,erp_data,0);
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
sig_report_fname = [SBJ_vars.dirs.stats SBJ '_mGLM_ROI_' model_id '_' stat_id '_' an_id '_sig_report.txt'];
if exist(sig_report_fname)
    system(['mv ' sig_report_fname ' ' sig_report_fname(1:end-4) '_bck.txt']);
end
sig_report = fopen(sig_report_fname,'a');
if st.rt_corr; reg_lab = [reg_lab {'RT'}]; end
result_str = ['%-10s' repmat('%-10i',[1 numel(reg_lab)]) '\n'];

% Print header
fprintf(sig_report,'%s (n = %i)\n',SBJ,numel(beta.label));
fprintf(sig_report,[repmat('%-10s',[1 1+numel(reg_lab)]) '\n'],'label',reg_lab{:});

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
sig_mat = zeros([numel(beta.label) numel(reg_lab)]);
for ch_ix = 1:numel(beta.label)
    % Consolidate to binary sig/non-sig
    for reg_ix = 1:numel(reg_lab)
        if any(squeeze(beta.qval(reg_ix,ch_ix,:))<st.alpha)
            sig_mat(ch_ix,reg_ix) = 1;
        end
    end
    if st.rt_corr && any(rt.mask(ch_ix,1,:))
        sig_mat(ch_ix,numel(reg_lab)+1) = 1;
    end
    
    % Report on significant electrodes for this SBJ
    fprintf(sig_report,result_str,beta.label{ch_ix},sig_mat(ch_ix,:));
end

fclose(sig_report);

%% Save Results
out_fname = [SBJ_vars.dirs.stats SBJ '_mGLM_ROI_' model_id '_' stat_id '_' an_id '.mat'];
if st.rt_corr
    save(out_fname,'-v7.3','beta','rt','st');
else
    save(out_fname,'-v7.3','beta','st');
end

end
