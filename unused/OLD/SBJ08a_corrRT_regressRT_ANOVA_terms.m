function SBJ08a_corrRT_regressRT_ANOVA_terms(SBJ,an_id,stat_id)
%% Run ANOVA with potential RT regression before
% Set up paths
% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

rng('shuffle'); % seed randi with time

%% Load Data
sbj_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(sbj_cmd);
an_cmd = ['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m'];
eval(an_cmd);
stat_cmd = ['run ' root_dir 'PRJ_Error/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_cmd);

load(strcat(SBJ_vars.dirs.events,SBJ,'_trl_info_final.mat'));
load(strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'.mat'));

% Select data in stat window
cfg_trim = [];
cfg_trim.latency = [stat_lim(1) stat_lim(2)+0.001]; %add a data point to get a full window over the tail end
hfa = ft_selectdata(cfg_trim,hfa);

%% Build Design Matrix
design = zeros([numel(trl_info.trl_n) numel(groups)]);
for grp_ix = 1:numel(groups)
    design(:,grp_ix) = fn_condition_index(groups{grp_ix}, trl_info)';
end

%% Run correlations with RT
if rt_correlation
    cfg_rt.design           = zscore(trl_info.rt);
    cfg_rt.ivar             = 1;
    stat = ft_freqstatistics(cfg_rt, hfa);
end
% OUTPUT:
%   .rho   = correlation coefficients
%   .stat  = t values
%   .prob  = p values

%% Regress off Reaction Time
if regress_rt
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
% NOTE: .stat and .prob output if .statistics = 'yes', but DO NOT TRUST THEM!!! (e.g., p values of 2...)

%% Run ANOVA
% Sliding window parameters
win_lim    = fn_sliding_window_lim(squeeze(hfa.powspctrm(1,1,1,:)),win_len,win_step);
win_center = round(mean(win_lim,2));

% Create structure for w2 in fieldtrip style
w2.cond   = {};
for factor_ix = 1:size(anova_terms,1)
    if sum(anova_terms(factor_ix,:))<2
        w2.cond = [w2.cond groups(anova_terms(factor_ix,:)~=0)];
    else
        w2.cond = [w2.cond strjoin(groups(anova_terms(factor_ix,:)~=0),'*')];
    end
end
w2.time   = hfa.time(win_center);
w2.label  = hfa.label;
w2.dimord = 'rpt_chan_time';
w2.trial  = zeros([numel(w2.cond) length(hfa.label) length(w2.time)]);
w2.boot   = zeros([numel(w2.cond) length(hfa.label) length(w2.time) n_boots]);
w2.pval   = w2.trial;

% Compute ANOVA and Explained Variance
for ch_ix = 1:length(hfa.label)
    fprintf('============================ %s (%i / %i) ============================\n',...
        hfa.label{ch_ix},ch_ix,numel(hfa.label));
    fprintf('   time = ');
    for t_ix = 1:numel(win_center);
        fprintf('%.3f, ',hfa.time(win_center(t_ix)));
        % Compute ANOVA (real 1st iteration, then bootstrap with randomized design matrix
        for boot_ix = 1:n_boots+1
            % Randomize design matrix after first real model
            grp_col = {};
            if boot_ix == 1 % Real Model
                for grp_ix = 1:numel(groups)
                    grp_col = {grp_col{:} design(:,grp_ix)};
                end
            else            % Randomized Model
                for grp_ix = 1:numel(groups)
                    grp_col = {grp_col{:} design(randi(size(design,1),size(design(:,grp_ix))),grp_ix)};
                end
            end
            [p, table] = anovan(squeeze(nanmean(hfa.powspctrm(:,ch_ix,1,win_lim(t_ix,1):win_lim(t_ix,2)),4)), ...
                grp_col, ...% design(:,2)}, ...
                'model', anova_terms, 'sstype', 2, ...% 'continuous', strmatch('RT',w2.cond),
                'varnames', groups, 'display', 'off');
            
            % Calculate w2 (debiased effect size; multiply with 100 to get PEV)
            %   table(:,2) = sum of squares
            %   table(:,3) = dof
            %   table(:,5) = mean square error
            %   table(:,6) = F statistic for main effects, only for rows = groups*2 (2nd set are interactions)
            %   table(:,7) = p value for main effects, only for rows = groups*2 (2nd set are interactions)
            %   rows: 1 = labels, 2:n_cond+2 = groups, n_cond+3 = DifOut interaction, end-1 = error, end = total
            %       interactions are ordered: 1*2, 1*3, ..., 1*n, 2*3, ..., 2*n, 3*4, ..., etc.
            %   w2 = (SSbw - df*MSE) / (SStot + MSE)
            mse = table{numel(w2.cond)+2,5};
            for factor_ix = 1:numel(w2.cond)
                factor_row = strmatch(w2.cond{factor_ix},table(:,1),'exact');
                if boot_ix == 1
                    w2.trial(factor_ix,ch_ix,t_ix) = (table{factor_row,2} - (table{factor_row,3} * mse))/...
                        (table{end,2} + mse);
                else
                    w2.boot(factor_ix,ch_ix,t_ix,boot_ix-1) = (table{factor_row,2} - (table{factor_row,3} * mse))/...
                        (table{end,2} + mse);
                end
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

%% Save Results
out_fname = [SBJ_vars.dirs.proc SBJ '_ROI_' stat_id '_' an_id '.mat'];
if rt_correlation
    save(out_fname,'-v7.3','hfa','w2','stat');
else
    save(out_fname,'-v7.3','hfa','w2');
end

end
