function SBJ08g_HFA_crRT_mLME(SBJs, proc_id, an_id, model_id, stat_id, atlas_id, roi_id, lme_formula)
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
    
    % select regressors (provisional)
%     model = model(:,1:end-1);
%     reg_lab = reg_lab(1:end-1);
%     disp(reg_lab)
    
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
    fprintf('==================== Converting to long format ========================\n');
    ntimes = size(hfa_win,3);
    
    for r = 1:numel(roi_list)
        croi = roi_list{r};
        if s == 1
            hfa_data{r} = cell(numel(SBJs), ntimes);
        end
        cur_chan = chan_labs(strcmp(roi_labs,croi));
        nch = numel(cur_chan);
        for ct = 1:ntimes
            hfa_data{r}{s,ct} = cell(size(hfa_win,1)*nch,3 + size(model,2));
            dcount = 0;
            for ch = 1:nch
                lab = cur_chan{ch};
                cdata = squeeze(hfa_win(:,ch,ct));
                for cdp = 1:numel(cdata)
                    dcount = dcount + 1;
                    hfa_data{r}{s,ct}{dcount,1} = cdata(cdp);
                    hfa_data{r}{s,ct}{dcount,2} = SBJ;
                    hfa_data{r}{s,ct}{dcount,3} = [SBJ '_' lab];
                    for mm = 1:size(model,2)
                        hfa_data{r}{s,ct}{dcount,3 + mm} = model(cdp,mm);
                    end
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
        ctable.Properties.VariableNames = [{'y','sub','chan'},reg_lab];
        hfa_tables{r}{ct} = ctable;
    end
end
clear hfa_data

%% Run LMEs
fprintf('==================== Fitting LME ========================\n');
LMEs = {};
coefs = {};
lower = {};
upper = {};
pvals = {};
rng(090590)
for r = 1:numel(hfa_tables)
    nt = size(hfa_tables{r},1);
    LMEs{r} = cell(nt,1);
    pvals{r} = NaN(size(model,2) + 1,nt);
    coefs{r} = NaN(size(model,2) + 1,nt);
    lower{r} = NaN(size(model,2) + 1,nt);
    upper{r} = NaN(size(model,2) + 1,nt);
    for t = 1:nt
        fprintf('fitting label %d and time point %d\n',r,t)
        
        ctable = hfa_tables{r}{t};
        LMEs{r}{t} = fitlme(ctable, lme_formula);%, 'FitMethod','REML');
        pvals{r}(:,t) = LMEs{r}{t}.Coefficients.pValue;
        coefs{r}(:,t) = LMEs{r}{t}.Coefficients.Estimate;
        lower{r}(:,t) = LMEs{r}{t}.Coefficients.Lower;
        upper{r}(:,t) = LMEs{r}{t}.Coefficients.Upper;
    end
end
%% FDR correction
fprintf('==================== FDR correction ========================\n');
qvals = cell(size(pvals));
for r = 1:numel(pvals)
    qvals{r} = NaN(size(pvals{r}));
    for cf = 1:size(pvals{r},1)
        [~, ~, ~, qvals{r}(cf,:)] = fdr_bh(squeeze(pvals{r}(cf,:)));
    end
end
% for r = 1:numel(qvals)
%     qvals{r} = NaN(size(pvals{r}));
%     qvals{r}(1,:) = pvals{r}(1,:); % don't include intercept in correction
%     [~, ~, ~, qvals{r}(2:end,:)] = fdr_bh(pvals{r}(2:end,:));
% end
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
beta.dimord = 'coef_time';

%% Save
fprintf('==================== Saving results ========================\n');
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
out_fname = [stats_dir model_id '_' stat_id '_' an_id '_hfa.mat'];
save(out_fname,'-v7.3','beta')

%% Extract electrode-wise coefficients.
chan_coef = {};
chan_lower = {};
chan_upper = {};
chan_pval = {};
chan_labels = {};
for r = 1:numel(LMEs)
    chan_coef{r} = NaN(size(model,2) + 1, length(unique(hfa_tables{r}{1}.chan)), numel(LMEs{r}));
    chan_lower{r} = chan_coef{r};
    chan_upper{r} = chan_coef{r};
    chan_pval{r} = chan_coef{r};
    for t = 1:numel(LMEs{r})
        [~,~,rndstats] =  randomEffects(LMEs{r}{t});
        
        %select channel random effect
        rndstats = rndstats(strcmp(rndstats.Group,'sub:chan'),:);
        if t == 1
            chan_labels{r} = unique(rndstats.Level);
        end
        chan_coef{r}(:,:,t) = reshape(rndstats.Estimate,size(chan_coef{r},1:2));
        chan_lower{r}(:,:,t) = reshape(rndstats.Lower,size(chan_coef{r},1:2));
        chan_upper{r}(:,:,t) = reshape(rndstats.Upper,size(chan_coef{r},1:2));
        chan_pval{r}(:,:,t) = reshape(rndstats.pValue,size(chan_coef{r},1:2));
    end
    chan_coef{r} = permute(chan_coef{r},[2,1,3]);
    chan_lower{r} = permute(chan_lower{r},[2,1,3]);
    chan_upper{r} = permute(chan_upper{r},[2,1,3]);
    chan_pval{r} = permute(chan_pval{r},[2,1,3]);
end

% FDR correction
chan_qval = {};
for r = 1:numel(chan_pval)
    for ch = 1:size(chan_pval{r},1)
        for rg = 1:size(chan_pval{r},2)
            [~,~,~,chan_qval{r}(ch,rg,:)] = fdr_bh(chan_pval{r}(ch,rg,:));
        end
    end
end

% Store in a structure
beta_chan =  [];
%conn_stats.LMEs = LMEs;
beta_chan.feature   = reg_lab;
beta_chan.time      = st.win_center;
beta_chan.coefs = chan_coef;
beta_chan.lower = chan_lower;
beta_chan.upper = chan_upper;
beta_chan.pvals = chan_pval;
beta_chan.qvals = chan_qval;
beta_chan.label = roi_list;
beta_chan.chan_label = chan_labels;
beta_chan.win_lim   = win_lim;
beta_chan.win_lim_s = hfa_time(win_lim);
beta_chan.dimord = 'chan_coef_time';

fprintf('=================== Saving channel effects ======================\n');
stats_dir = [root_dir 'PRJ_Error/data/GRP/stats/'];
out_fname = [stats_dir model_id '_' stat_id '_' an_id '_hfa_chancoef.mat'];
save(out_fname,'-v7.3','beta_chan')
