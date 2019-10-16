function SBJ08a_crRT_nANOVA(SBJ,an_id,stat_id)
%% function SBJ08a_crRT_mANOVA(SBJ,an_id,stat_id)
%   Run ANOVA with sliding windows for given time-frequency analysis
%   Correlation with RT and regression of RT as confound supported
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
eval(['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m']);

load(strcat(SBJ_vars.dirs.events,SBJ,'_trl_info_final.mat'));
load(strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'.mat'));

% Check if more than one frequency, error for now
if numel(hfa.freq)>1
    error('HFA has more than one frequency, can''t run on that for now!');
end

% Select data in stat window
cfg_trim = [];
cfg_trim.latency = [st.stat_lim(1) st.stat_lim(2)+0.001];
hfa = ft_selectdata(cfg_trim,hfa);

%% Build Design Matrix
[grp_lab, ~, ~] = fn_group_label_styles(st.model_lab);
design = cell([1 numel(st.groups)]);
levels = cell([1 numel(st.groups)]);
for grp_ix = 1:numel(st.groups)
    [levels{grp_ix}, ~, ~] = fn_condition_label_styles(st.groups{grp_ix});
    design{grp_ix} = fn_condition_index(st.groups{grp_ix}, trl_info);
end

%% Exclude trials based on RT and trial type
% % Select trial types of interest
% if ~strcmp(st.trial_cond,'all')
%     good_cond_idx = zeros(size(trl_info.trial_n));
%     for cond_ix = 1:numel(st.trial_cond)
%         cond_idx = fn_condition_index(st.trial_cond{cond_ix}, trl_info.condition_n, 'trl_info', trl_info);
%         good_cond_idx(cond_idx) = 1;
%     end
% else
%     good_cond_idx = true(size(trl_info.trial_n));
% end
% 
% % Log combined bad trial types
% good_trl_idx = all([~bad_rt_idx good_cond_idx],2);
% st.bad_trials.all  = trl_info.trial_n(~good_trl_idx);
% st.bad_trials.rt   = trl_info.trial_n(bad_rt_idx);
% st.bad_trials.cond = trl_info.trial_n(~good_cond_idx);
% 
% % Exclude bad trials
% for grp_ix = 1:numel(st.groups)
%     design{grp_ix} = design{grp_ix}(good_trl_idx);
% end
% ti_fields = fieldnames(trl_info);
% orig_n_trials = numel(trl_info.trial_n);
% for f_ix = 1:numel(ti_fields)
%     if numel(trl_info.(ti_fields{f_ix}))==orig_n_trials
%         trl_info.(ti_fields{f_ix}) = trl_info.(ti_fields{f_ix})(good_trl_idx);
%     end
% end

%% Compute max z score per condition (for later exclusions)
max_z = zeros(size(hfa.label));
cfg_avg = [];
cfg_avg.avgoverrpt = 'yes';
for grp_ix = 1:numel(st.groups)
    for level_ix = 1:numel(levels{grp_ix})
        % Average within condition
        cfg_avg.trials = find(design{grp_ix}==level_ix);
        tmp_avg = ft_selectdata(cfg_avg,hfa);
        % Binarize channels by average HFA vs. threshold
        max_z = max(max_z,max(squeeze(tmp_avg.powspctrm),[],2));
    end
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
w2.design    = design;
w2.cond      = grp_lab;
w2.time      = hfa.time(win_center);
w2.win_lim   = win_lim;
w2.label     = hfa.label;
w2.dimord    = 'rpt_chan_time';
w2.max_hfa_z = max_z;
% w2.trial     = zeros([numel(w2.cond) length(hfa.label) length(w2.time)]);
w2.boot      = zeros([numel(w2.cond) length(hfa.label) length(w2.time) st.n_boots]);
% w2.pval      = w2.trial;
w2.win_lim_s = hfa.time(win_lim);

% Compute ANOVA and Explained Variance
for ch_ix = 1:numel(w2.label)
    fprintf('============================ %s (%i / %i) ============================\n',...
        w2.label{ch_ix},ch_ix,numel(w2.label));
    fprintf('   time = ');
    for t_ix = 1:numel(w2.time)
        fprintf('%.3f, ',w2.time(t_ix));
        % Compute ANOVA (real 1st iteration, then bootstrap with randomized design matrix
        [~, table] = anovan(squeeze(hfa_win(:,ch_ix,t_ix)), w2.design, ...
            'model', st.anova_term,...% 'sstype', 2, ...% 'continuous', strmatch('RT',w2.cond),
            'varnames', st.groups, 'display', 'off');
        
        % Calculate w2 (debiased effect size; multiply with 100 to get PEV)
        %   rows: 1 = labels, 2:n_cond+2 = groups, n_cond+3 = DifOut interaction, end-1 = error, end = total
        %       interactions are ordered: 1*2, 1*3, ..., 1*n, 2*3, ..., 2*n, 3*4, ..., etc.
        %   w2 = (SSbw - df*MSE) / (SStot + MSE)
        sumsq_ix = find(strcmp(table(1,:),'Sum Sq.'));
        dof_ix   = find(strcmp(table(1,:),'d.f.'));
        mse_ix   = find(strcmp(table(1,:),'Mean Sq.'));
%         fstat_ix = find(strcmp(table(1,:),'F'));
%         pval_ix  = find(strcmp(table(1,:),'Prob>F'));
        mse = table{numel(w2.cond)+2,mse_ix};
        for factor_ix = 1:numel(w2.cond)
            factor_row = strcmp(w2.cond{factor_ix},table(:,1));
            w2.trial(factor_ix,ch_ix,t_ix) = (table{factor_row,sumsq_ix} - (table{factor_row,dof_ix} * mse))/...
                (table{end,sumsq_ix} + mse);
        end
        clear p table
        
        % Compute ANOVA for permuted data
        rand_design = design;
        fprintf('boot #: ');
        for boot_ix = 1:st.n_boots
            fprintf('%i..',boot_ix);
            % Randomize rows in design matrix
            rand_idx = randperm(size(design{grp_ix},1));
            for grp_ix = 1:numel(st.groups)
                rand_design{grp_ix} = design{grp_ix}(rand_idx);
            end
            % Compute ANOVA
            [~, table] = anovan(squeeze(hfa_win(:,ch_ix,t_ix)), w2.design, ...
                'model', st.anova_term,...% 'sstype', 2, ...% 'continuous', strmatch('RT',w2.cond),
                'varnames', st.groups, 'display', 'off');
            % Compute w2
            mse = table{numel(w2.cond)+2,mse_ix};
            for factor_ix = 1:numel(w2.cond)
                factor_row = strcmp(w2.cond{factor_ix},table(:,1));
                w2.boot(factor_ix,ch_ix,t_ix,boot_ix) = (table{factor_row,sumsq_ix} - (table{factor_row,dof_ix} * mse))/...
                    (table{end,sumsq_ix} + mse);
            end
            if mod(boot_ix,100)==0
                fprintf('\n');
            end
            clear p table
        end
        
        % Compute Significance accross bootstraps
        for factor_ix = 1:numel(w2.cond)
            w2_false_pos = find(w2.boot(factor_ix,ch_ix,t_ix,:) > w2.trial(factor_ix,ch_ix,t_ix));
            w2.pval(factor_ix,ch_ix,t_ix) = numel(w2_false_pos)/n_boots;
        end
    end
    fprintf('\n');
end

% Compute statistics
% w2.pval = sum(bsxfun(@ge,w2.boot,w2.trial),4)./st.n_boots; % sum(boots>real)/n_boots
w2.zscore   = norminv(1-w2.pval,0,1);
w2.bootmean = mean(w2.boot,4);
w2.bootstd  = std(w2.boot,[],4);
% w2 = rmfield(w2,'boot');
w2.zscore(isinf(w2.zscore)) = norminv(1-1/st.n_boots/2,0,1);

% Multiple Comparisons Correction within Channel
w2.qval = nan(size(w2.pval));
for ch_ix = 1:numel(w2.label)
    [~, ~, ~, w2.qval(:,ch_ix,:)] = fdr_bh(squeeze(w2.pval(:,ch_ix,:)));
end

%% Print results
% Prep report
sig_report_fname = [SBJ_vars.dirs.proc SBJ '_nANOVA_ROI_' stat_id '_' an_id '_sig_report.txt'];
if exist(sig_report_fname)
    system(['mv ' sig_report_fname ' ' sig_report_fname(1:end-4) '_bck.txt']);
end
sig_report = fopen(sig_report_fname,'a');
if st.rt_corr; cond_lab = [grp_lab {'RT'}]; else cond_lab = grp_lab; end
result_str = ['%-10s' repmat('%-10i',[1 numel(cond_lab)]) '\n'];

% Print header
fprintf(sig_report,'%s (n = %i)\n',SBJ,numel(w2.label));
fprintf(sig_report,[repmat('%-10s',[1 1+numel(cond_lab)]) '\n'],'label',cond_lab{:});

% Print summary lines (absolute)
if st.rt_corr
    fprintf(sig_report,result_str, 'count', sum(any(w2.qval(:,:,:)<st.alpha,3),2), sum(any(rt.mask(:,1,:),3)));
    fprintf(sig_report,strrep(result_str,'i','.3f'), 'percent',...
        [sum(any(w2.qval(:,:,:)<st.alpha,3),2)', sum(any(rt.mask(:,1,:),3))]./numel(w2.label));
else
    fprintf(sig_report,result_str, 'count', sum(any(w2.qval(:,:,:)<st.alpha,3),2));
    fprintf(sig_report,strrep(result_str,'i','.3f'), 'percent',...
        sum(any(w2.qval(:,:,:)<st.alpha,3),2)'./numel(w2.label));
end

% Print Channel Lines
sig_mat = zeros([numel(w2.label) numel(cond_lab)]);
for ch_ix = 1:numel(w2.label)
    % Consolidate to binary sig/non-sig
    for grp_ix = 1:numel(grp_lab)
        if any(squeeze(w2.qval(grp_ix,ch_ix,:))<st.alpha)
            sig_mat(ch_ix,grp_ix) = 1;
        end
    end
    if st.rt_corr && any(rt.mask(ch_ix,1,:))
        sig_mat(ch_ix,grp_ix) = 1;
    end
    
    % Report on significant electrodes for this SBJ
    fprintf(sig_report,result_str,w2.label{ch_ix},sig_mat(ch_ix,:));
end

fclose(sig_report);

%% Save Results
out_fname = [SBJ_vars.dirs.proc SBJ '_nANOVA_ROI_' stat_id '_' an_id '.mat'];
if st.rt_corr
    save(out_fname,'-v7.3','w2','rt','st');
else
    save(out_fname,'-v7.3','w2','st');
end

end
