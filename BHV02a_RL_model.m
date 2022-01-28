function BHV02a_RL_model(SBJ,proc_id,model_id)
%% Run reinforcement learning (RL) model on single SBJ behavior
%   for specific combination of model regressors and conditions
% Plots and saves:
%   design matrix (SBJ_modelID_design.png)
%   correlation matrix between regressors (SBJ_modelID_design_corr.png)
%   violin plots of regressors by condition (SBJ_modelID_reg_cond.png)
% Expected value computation excludes neutral (surprise) outcomes
% Accounts for paradigm restarts (e.g., after crash), excludes training
% -------------------------------------------------------------------------
% WARNING: bhv.sound was added during revisions, and SBJ02a and SBJ02c
%   were not re-run to propagate that field to SBJ_behav_final.mat.
%   Instead, trial sounds are taken directly via fn_load_behav_csv here
%   when loading the original trial list, then added to the model in the
%   right order.
% INPUTS:
%   SBJ [str] - ID of subject to run
%   proc_id [str] - ID of preprocessing pipeline
%   model_id [str] - ID of the model parameters to use
%       mdl.model_lab [str] - selects set of regressors for model
%       mdl.model_cond [cell string] - selects conditions to model
%       mdl.bias_reg [cell array] - str names of regressor to bias
%       mdl.bias [array] - amounts to add to biased regressor
%       mdl.bias_cond [cell array] - strs selecting which conditions to bias
% OUTPUTS:
%   model [float array] - [n_trials, n_regressors] matrix of all regressors
%   ========== Main RL Model Regressors ==========
%   Lik: likelihood/probability of given outcome in given difficulty context
%       rLik: same but residualized (against uRPE)
%   EV: expected value, linearly scaled from probability to -1:1 reward range
%       pWin: win probabilty from logistic regression
%       bAcc: block level accuracy
%       rAcc: rolling accuracy (last 5 trials; can be changed via roll_win;
%           accuracy computed within difficulty context: roll_wi_cond = true)
%   sRPE: signed reward prediction error (RPE)
%   uRPE: unsigned reward prediction error (RPE) = abs(sRPE)
%   ========== Unused Outcome Regressors ==========
%   Val: reward value (-1/0/1)
%   Mag: reward magnitude = abs(Val)
%   Sign: valence according to RPE (accounts for expectations, so neutral
%       is negative in easy and positive in hard)
%   ========== Extra Unused Control Regressors ==========
%   score: cumulative score across the experiment
%   ITI: inter-trial interval on the previous trial
%   ========== Extra Auditory Salience Regressors ==========
%   ERB: "equivalent rectangular bandwidth", roughly loudness
%   roughness: Sijia 2019 JNeuro- strongly correlates with subjective
%       salience and physiological index of saleince (pupil dilation)
%   maxSalience: max amplitude in Kayser (2006) salience map model
%   meanSalience: mean amplitude in Kayser (2006) salience map model
%   maxGradient:  max change in time in Kayser (2006) salience map model
%   meanGradient: mean change in time in Kayser (2006) salience map model
%   ========== Unused Performance (RT) Regressors ==========
%       These are all log transformed since humans perceive distance in log scale
%       Precision decreases with distance, error = precision^-1
%   sTaEr: signed target error (distance from RT to target time)
%   sTaPr: signed target precision
%   uTaEr: unsigned target error
%   uTaPr: unsigned target precision
%   old options included previous trial sTaEr
%   sThEr: signed threshold error (distance from RT to nearest tolerance threshold)
%   sThPr: signed threshold precision
%   uThEr: unsigned threshold error
%   uThPr: unsigned threshold precision
% betas [float array] - coefficeints from logistic regression fit

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/hoycw/';app_dir=[root_dir 'Apps/'];
else; root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Load Data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
model_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/model_vars/' model_id '_vars.m'];
eval(model_vars_cmd);

% Determine model parameteres and conditions
if strcmp(mdl.model_lab,'RL3D')
    [reg_lab, reg_names, ~, ~] = fn_performanceRL_regressor_label_styles(mdl.model_lab);
else
    [reg_lab, reg_names, ~, ~, ~] = fn_regressor_label_styles(mdl.model_lab);
end
[cond_lab, cond_names, cond_colors, ~, ~] = fn_condition_label_styles(mdl.model_cond);

% Check bias parameters
if isfield(mdl,'bias') || isfield(mdl,'bias_reg') || isfield(mdl,'bias_cond')
    % Check all necessary fields exist
    if ~isfield(mdl,'bias') || ~isfield(mdl,'bias_reg') || ~isfield(mdl,'bias_cond')
        error([model_id ' is missing some bias fields!']);
    end
    % Check all necessary fields are same length
    if numel(mdl.bias_reg)~=numel(mdl.bias) || numel(mdl.bias_reg)~=numel(mdl.bias_cond)
        error([model_id ' has mismatched number of bias terms!']);
    end
    % Only allow bias for pWin for now
    if any(~strcmp(mdl.bias_reg,'pWin'))
        error('Bias is only set up for pWin right now!');
    end
end

% Check residualization parameters
if isfield(mdl,'resid')
    error('Decided to not residualize for now, so do not use this!');
    % Avoid order effects by limiting to one residalization
    if numel(mdl.resid)>1; error('Only allowing one variable to be residualized for now'); end
    % Check both regressors exist
    for r_ix = 1:numel(mdl.resid)
        if ~any(strcmp(reg_lab,mdl.resid{r_ix}{1})) || ~any(strcmp(reg_lab,mdl.resid{r_ix}{2}))
            error(['residualize regressor missing: ' strjoin(mdl.resid{r_ix},',')]);
        end
    end
end

%% Load and Select Behavior
% Load data
load([root_dir 'PRJ_Error/data/' SBJ '/03_events/' SBJ '_bhv_' proc_id '_final.mat']);

% Select Conditions of Interest
full_cond_idx = fn_condition_index(cond_lab, bhv);
full_n_trials = numel(bhv.trl_n);
bhv = fn_select_bhv(bhv, full_cond_idx);
n_trials = numel(bhv.trl_n);
fprintf('\t%s: Loaded %d trials, kept %d for modeling...\n',SBJ,full_n_trials,n_trials);

%% Load full trial history
% This section is unnecessary unless I bring back the rolling accuracy
%   control condion, which requires combining separate runs

% % Use all trials from behavior
% [bhv_orig] = fn_load_behav_csv([SBJ_vars.dirs.events SBJ '_behav.csv']);
% orig_n_trials = numel(bhv_orig.trl_n);
% 
% % Find fully visible training trials (at very beginning)
% full_vis_trl_n = 1:bhv.prdm.n_examples;
% full_vis_ix = [];
% for t_ix = 1:numel(full_vis_trl_n)
%     full_vis_ix = [full_vis_ix find(bhv_orig.trl_n==full_vis_trl_n(t_ix))];
% end
% full_vis_ix = sort(full_vis_ix(:));
% 
% % Find training trials
% train_ix = find(bhv_orig.blk==0);
% 
% % Remove full vis and training
% %   Keep training but not fully visible trials to compute rolling accuracy
% %   for first couple trials
% bhv_orig_nfv = bhv_orig;
% bhv_orig_ntr = bhv_orig;
% for f_ix = 1:numel(bhv_fields)
%     if numel(bhv_orig.(bhv_fields{f_ix}))==orig_n_trials
%         bhv_orig_nfv.(bhv_fields{f_ix}) = bhv_orig.(bhv_fields{f_ix})(setdiff(1:numel(bhv_orig.trl_n),full_vis_ix));
%         bhv_orig_ntr.(bhv_fields{f_ix}) = bhv_orig.(bhv_fields{f_ix})(setdiff(1:numel(bhv_orig.trl_n),train_ix));
%     end
% end

%% Compute Reward Value and Magnitude
% Reward value (Val) = score change for given outcome
%   Rewards are scaled from -100/0/100 to -1/0/1 for simplicity
% Reward magnitude (Mag) = absolute value of reward value
% Note this is without accounting for expected value
if any(strcmp(reg_lab,'Val')) || any(strcmp(reg_lab,'Mag'))
    Val = nan(size(bhv.trl_n));
    Mag = nan(size(bhv.trl_n));
    for cond_ix = 1:numel(cond_lab)
        idx = logical(fn_condition_index(cond_lab(cond_ix),bhv));
        % Reward value
        switch cond_lab{cond_ix}
            case 'EzWn'
                Val(idx) = 1;
            case 'EzLs'
                Val(idx) = -1;
            case 'EzNu'
                Val(idx) = 0;
            case 'HdWn'
                Val(idx) = 1;
            case 'HdLs'
                Val(idx) = -1;
            case 'HdNu'
                Val(idx) = 0;
            otherwise
                error(['Unknown trial type: ' cond_lab{cond_ix}]);
        end
        % Reward magnitude
        Mag(idx) = abs(Val(idx));
    end
end

%% Compute Outcome Likelihood
% Outcome likelihood (Lik) = proportion of trials for given outcome in
%   given difficulty context over the whole experiment (except training)
% This regressor is referred to as 'probability' in the paper
if any(strcmp(reg_lab,'Lik')) || any(strcmp(reg_lab,'rLik'))
    Lik = nan(size(bhv.trl_n));
    for cond_ix = 1:numel(cond_lab)
        idx = fn_condition_index(cond_lab(cond_ix),bhv);
        % Compute within difficulty context, which is known to participant
        if strcmp(cond_lab{cond_ix}(1:2),'Ez')
            n_trls = sum(strcmp('easy',bhv.cond));
        elseif strcmp(cond_lab{cond_ix}(1:2),'Hd')
            n_trls = sum(strcmp('hard',bhv.cond));
        else
            error('Neither Ez nor Hd, what is it?');
        end
        Lik(logical(idx)) = sum(idx)/n_trls;
    end
    if any(strcmp(reg_lab,'rLik'))
        rLik = Lik;
    end
end

%% Compute Win Probabilities and Expected Value
% EV: Logistic regression to predict win/loss based on tolerance
% bAcc: block-level accuracy
% rAcc: rolling accuracy in last 5 or 10 trials
% Probability of winning is linearly scaled to range of reward value
% Exclude neutral trials since feedback is not determined by performance
n_idx = fn_condition_index({'Nu'},bhv);
if any(strcmp(reg_lab,'EV')) || any(strcmp(reg_lab,'Sign')) ...
        || any(strcmp(reg_lab,'sRPE')) || any(strcmp(reg_lab,'uRPE'))
    % Select Data (fit on everything except neutral)
    X = bhv.tol(~n_idx);
    y = double(bhv.hit(~n_idx));
    
    % Logistic regression
    betas = glmfit(X,y,'binomial','link','logit');
    
    % Generated probability of winning from model
    z = betas(1) + (bhv.tol * betas(2));
    pWin = 1 ./ (1+exp(-z));
    
    % Add bias term
    if isfield(mdl,'bias_reg') && any(strcmp(mdl.bias_reg,'pWin'))
        for b_ix = find(strcmp(mdl.bias_reg,'pWin'))
            bias_idx = fn_condition_index(mdl.bias_cond(b_ix),bhv);
            pWin(bias_idx~=0) = pWin(bias_idx~=0) + mdl.bias(b_ix);
        end
    end
    
    % Linearly scale from probability (0 to 1) to expected value (-1 to 1)
    EV = pWin*2 - 1;
end

% Alternative computations of expected value
if any(strcmp(reg_lab,'bAcc'))
    % Get block level accuracy per run/block
    blk_acc = nan([numel(unique(bhv.run)) numel(unique(bhv.blk))]);
    for r_ix = 1:numel(unique(bhv.run))
        run_blks = unique(bhv.blk(bhv.run==r_ix));
        for b_ix = 1:numel(run_blks)
            blk_acc(r_ix,run_blks(b_ix)) = mean(strcmp(bhv.fb(bhv.run==r_ix & bhv.blk==run_blks(b_ix) & ~n_idx),'W'));
        end
    end
    
    % Format block accuracy as single trial predictor
    bAcc = nan(size(bhv.trl_n));
    for t_ix = 1:numel(bhv.trl_n)
        bAcc(t_ix) = blk_acc(bhv.run(t_ix), bhv.blk(t_ix));
    end
    
    % Linearly scale from probability (0 to 1) to expected value (-1 to 1)
    EV = bAcc*2 - 1;
elseif any(strcmp(reg_lab,'rAcc'))
    error('need to concatenate original bhv files for this');
%     % Analysis Parameters (to be finalized)
%     roll_win = 5;           % number of previous trials to average accuracy over
%     roll_wi_cond = true;    % only average over trials within same difficulty context
%     
%     % Build index of feedback
%     % Includes non-fully visible training to estimate for first few trials
%     fb_idx = nan(size(bhv_orig_nfv.trl_n));
%     fb_idx(strcmp(bhv_orig_nfv.fb,'W')) = 1;
%     fb_idx(strcmp(bhv_orig_nfv.fb,'L')) = 0;
%     
%     % Compute rolling accuracy of last 5 trials
%     %   Include training trials!
%     run_idx_nfv = ones(size(bhv_orig_nfv.trl_n));
%     if any(strcmp(SBJ,{'EP09','EEG12','EEG24','EEG25'}))
%         new_run_ix_nfv = find(abs(diff(bhv_orig_nfv.trl_n))>bhv.prdm.n_trials)+1;
%         run_idx_nfv(new_run_ix_nfv:end) = 2;
%     end
%     
%     rAcc = nan(size(bhv.trl_n));
%     for t_ix = 1:numel(bhv.trl_n)
%         if roll_wi_cond
%             prior_cond_ix = find(strcmp(bhv_orig_nfv.cond,bhv.cond(t_ix)) & bhv_orig_nfv.trl_n<bhv.trl_n(t_ix));
%             prior_cond_ix = prior_cond_ix(run_idx_nfv(prior_cond_ix)==run_idx(t_ix));
%             rAcc(t_ix) = nanmean(fb_idx(prior_cond_ix(end-roll_win+1:end)));
%         else
%             error('not ready to compute rolling accuracy across conditions in double runs for EEG12,24,25');
%             orig_ix = find(bhv_orig_nfv.trl_n==bhv.trl_n(t_ix));
%             rAcc(t_ix) = nanmean(fb_idx(orig_ix-roll_win:orig_ix-1));
%         end
%     end
%     
%     % Linearly scale from probability (0 to 1) to expected value (-1 to 1)
%     EV = rAcc*2 - 1;
end

%% Compute Reward Valence (sign)
% valence of RPE (accounts for expectations)
if any(strcmp(reg_lab,'Sign'))
    Sign = sign(double(bhv.score)/100 - EV);
end

%% Compute Prediction Errors
% RPE = actual reward value - expected value
if any(strcmp(reg_lab,'sRPE'))
    sRPE = double(bhv.score)/100 - EV;
end
if any(strcmp(reg_lab,'uRPE'))
    uRPE = abs(double(bhv.score)/100 - EV);
end

if any(strcmp(reg_lab,'pRPE'))
    sRPE = double(bhv.score)/100 - EV;
    pRPE = zeros(size(sRPE));
    pRPE(sign(sRPE)==1) = sRPE(sign(sRPE)==1);
end
if any(strcmp(reg_lab,'nRPE'))
    sRPE = double(bhv.score)/100 - EV;
    nRPE = zeros(size(sRPE));
    nRPE(sign(sRPE)==-1) = -1*sRPE(sign(sRPE)==-1); % flip sign so positive coefficient means increasing activity for more negative RPE
end

%% Compute total, cumulative score
if any(strcmp(reg_lab,'score'))  
    error('need to concatenate original bhv files for this');
%     % Compute total score (exclude training trials)
%     full_score = zeros(size(bhv_orig_ntr.score));
%     new_run_ix_ntr = find(abs(diff(bhv_orig_ntr.trl_n))>bhv.prdm.n_trials)+1;
%     first_trl_ix = 1;
%     for t_ix = 1:numel(bhv_orig_ntr.trl_n)
%         % If paradigm restarted, sum over new session score (as presented on screen)
%         if n_runs==2 && t_ix==new_run_ix_ntr
%             first_trl_ix = new_run_ix_ntr;
%         end
%         full_score(t_ix) = sum(bhv_orig_ntr.score(first_trl_ix:t_ix));
%     end
%     
%     % Select trials
%     score = zeros(size(bhv.trl_n));
%     run_ix = 1;
%     for t_ix = 1:numel(bhv.trl_n)
%         if n_runs==2 && t_ix==new_run_ix
%             run_ix = 2;
%         end
%         ntr_ix = find(bhv_orig_ntr.trl_n==bhv.trl_n(t_ix));
%         if numel(ntr_ix)>1
%             ntr_ix = ntr_ix(run_ix);
%         end
%         score(t_ix) = full_score(ntr_ix);
%     end
end

%% Compute Auditory Salience Features (via Sijia Zhao)
% Check if any auditory salience regressors are requested
if ~isempty(intersect(reg_lab,{'ERB','rough','mxS','mnS','mxdS','mndS'}))
    % Load auditory salience features
    aud_sal = fn_load_auditory_salience([root_dir 'PRJ_Error/scripts/auditory_salience_output.csv']);
    
    % Initialize variables
    sal_reg_lab = intersect(reg_lab,{'ERB','rough','mxS','mnS','mxdS','mndS'});
    for reg_ix = 1:numel(sal_reg_lab)
        eval([sal_reg_lab{reg_ix} ' = nan(size(bhv.trl_n));']);
    end
    
    for trl_ix = 1:numel(bhv.trl_n)
        if strcmp(bhv.fb{trl_ix},'N')
            sound_name = [bhv.sound{trl_ix}(strfind(bhv.sound{trl_ix},'/')+1:end) '.wav'];
        else
            sound_name = [bhv.sound{trl_ix} '.wav'];
        end
        
        % ERB "loudness"
        if any(strcmp(reg_lab,'ERB'))
            ERB(trl_ix) = aud_sal.ERB(strcmp(aud_sal.sound,sound_name));
        end
        
        % Roughness
        if any(strcmp(reg_lab,'rough'))
            rough(trl_ix) = aud_sal.rough(strcmp(aud_sal.sound,sound_name));
        end
        
        % Max Salience
        if any(strcmp(reg_lab,'mxS'))
            mxS(trl_ix) = aud_sal.mxS(strcmp(aud_sal.sound,sound_name));
        end
        
        % Mean Salience
        if any(strcmp(reg_lab,'mnS'))
            mnS(trl_ix) = aud_sal.mnS(strcmp(aud_sal.sound,sound_name));
        end
        
        % Max Gradient (delta salience)
        if any(strcmp(reg_lab,'mxdS'))
            mxdS(trl_ix) = aud_sal.mxdS(strcmp(aud_sal.sound,sound_name));
        end
        
        % Mean Gradient (delta salience)
        if any(strcmp(reg_lab,'mndS'))
            mndS(trl_ix) = aud_sal.mndS(strcmp(aud_sal.sound,sound_name));
        end
    end
end

%% Distance from Target
% Signed Target Distance (RT - target_time to center for early/late)
% Unsigned version is just total error without early/late information
% Log transformed distance between RT and target time to account for humans
%   perceiving distances in log space
% Target Error (TaEr):  larger values for larger distances
% Target Precision (TaPr): smaller values for larger distances
tar = bhv.rt-bhv.prdm.target;
if any(strcmp(reg_lab,'sTaEr'))
    % Separately for wins and losses to maintain correct sign
    sTaEr = nan(size(tar));
    sTaEr(tar<0)  = log(-tar(tar<0));
    sTaEr(tar>=0) = -log(tar(tar>=0));
    % Invert to make large distances biggest values
    sTaEr = sTaEr.^-1;
end
if any(strcmp(reg_lab,'sTaPr'))
    % Separately for wins and losses to maintain correct sign
    sTaPr = nan(size(tar));
    sTaPr(tar<0)  = log(-tar(tar<0));
    sTaPr(tar>=0) = -log(tar(tar>=0));
end

% Unsigned Target Distance
if any(strcmp(reg_lab,'uTaEr'))
    % Flip sign to make it positive
    uTaEr = -log(abs(tar));
    % Invert to make large distances biggest values
    uTaEr = uTaEr.^-1;
end
if any(strcmp(reg_lab,'uTaPr'))
    % Flip sign to make it positive
    uTaPr = -log(abs(tar));
end

% Previous trial RTs
if any(strcmp(reg_lab,'psTaEr'))
    error('not using previous trial right now!');
    psTaEr = nan(size(bhv.rt));
    for t_ix = 2:numel(bhv.rt)
        psTaEr(t_ix) = sTaEr(t_ix-1);
    end
end
if any(strcmp(reg_lab,'p2sTaEr'))
    error('not using previous trial right now!');
    p2sTaEr = nan(size(bhv.rt));
    for t_ix = 3:numel(bhv.rt)
        p2sTaEr(t_ix) = psTaEr(t_ix-2);
    end
end

%% Distance from Tolerance Threshold
% Compute distance from closest tolerance bound (early vs. late RTs)
%   Positive = win, Negative = loss
% Can also be error or precision
% NOTE: Due to a rounding error in python paradigm code, very small percentage
%   of trials received the wrong feedback; since difference is less than
%   10ms, this error is likely imperceptable to SBJ, but those trials are
%   excluded from signed (but not unsigned) RT-threshold predictors
early_idx = (bhv.rt-bhv.prdm.target)<0;
thr = nan(size(bhv.rt));
thr(early_idx) = bhv.rt(early_idx)-(bhv.prdm.target-bhv.tol(early_idx));
thr(~early_idx) = (bhv.prdm.target+bhv.tol(~early_idx))-bhv.rt(~early_idx);

% Transform almost hit/miss into large values
if any(strcmp(reg_lab,'sThEr'))
    % Separately for wins and losses to maintain correct sign for beta
    % interpretation (+ for win, - for loss)
    sThEr = nan(size(bhv.rt));
    sThEr(bhv.hit==1) = -log(thr(bhv.hit==1));
    sThEr(bhv.hit==0) = log(-thr(bhv.hit==0));
    sThEr = sThEr.^-1;
    % Remove bad feedback trials in sThPr
    %   (incorrect feedback given due to rounding error)
    sThEr(logical(bhv.bad_fb)) = nan;
end
if any(strcmp(reg_lab,'sThPr'))
    % Separately for wins and losses to maintain correct sign for beta
    % interpretation (+ for win, - for loss)
    sThPr = nan(size(bhv.rt));
    sThPr(bhv.hit==1) = -log(thr(bhv.hit==1));
    sThPr(bhv.hit==0) = log(-thr(bhv.hit==0));
    % Remove bad feedback trials in sThPr
    %   (incorrect feedback given due to rounding error)
    sThPr(logical(bhv.bad_fb)) = nan;
end

% Add unsigned version
if any(strcmp(reg_lab,'uThEr'))
    % Separately for wins and losses to maintain correct sign for beta
    % interpretation (+ for win, - for loss)
    uThEr = nan(size(bhv.rt));
    uThEr(bhv.hit==1) = -log(thr(bhv.hit==1));
    uThEr(bhv.hit==0) = log(-thr(bhv.hit==0));
    uThEr = uThEr.^-1;
    % Keeping bad feedback trials because absolute distance is correct
    uThEr = abs(uThEr);
end
if any(strcmp(reg_lab,'uThPr'))
    % Separately for wins and losses to maintain correct sign for beta
    % interpretation (+ for win, - for loss)
    uThPr = nan(size(bhv.rt));
    uThPr(bhv.hit==1) = -log(thr(bhv.hit==1));
    uThPr(bhv.hit==0) = log(-thr(bhv.hit==0));
    % Keeping bad feedback trials because absolute distance is correct
    uThPr = abs(uThPr);
end

%% Inter-Trial Interval (from previous trial
ITI = bhv.ITI_type;
% Remove ITI on first trial of the block due to extended delay
ITI(ITI==0) = nan;
if any(~isnan(ITI(bhv.blk_trl_n==1)))
    error('First trial in block is not NaN in ITI regressor!');
end

%% Load Data and Build Model
% Concatenate all regressors into one design matrix
model = nan([sum(n_trials) numel(reg_lab)]);
for reg_ix = 1:numel(reg_lab)
    model(:,reg_ix) = eval(reg_lab{reg_ix});
end
% Report missing values
if any(isnan(model(:)))
    fprintf(2,'\tWARNING: %d NaNs detected in model!\n',sum(isnan(model(:))));
end

%% Residualize regressors
if isfield(mdl,'resid')
    for r_ix = 1:numel(mdl.resid)
        x_ix = find(strcmp(reg_lab,mdl.resid{r_ix}{1}));
        y_ix = find(strcmp(reg_lab,mdl.resid{r_ix}{2}));
        betas = polyfit(model(:,x_ix),model(:,y_ix),1);
        y_fit = polyval(betas,model(:,x_ix));
        model(:,x_ix) = model(:,y_ix) - y_fit;
    end
end

%% Compute and plot correlations between regressors
reg_corr = corr(model,'rows','complete');

% Create figure directory
fig_dir = [SBJ_vars.dirs.models mdl.model_id '_plots/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Plot design matrix
fig_name = [SBJ '_' mdl.model_id '_design'];
figure('Name',fig_name);
imagesc(model);
xticklabels(reg_lab);
colorbar;
saveas(gcf,[fig_dir fig_name '.png']);

% Plot regressor correlation matrix
fig_name = [SBJ '_' mdl.model_id '_design_corr'];
figure('Name',fig_name);
imagesc(reg_corr);
xticklabels(reg_lab);
yticklabels(reg_lab);
colorbar;
saveas(gcf,[fig_dir fig_name '.png']);

%% Plot Regressors by Condition
cond_idx = fn_condition_index(cond_lab, bhv);

fig_name = [SBJ '_' mdl.model_id '_reg_cond'];
figure('Name',fig_name,'units','normalized','outerposition',[0 0 1 1]);
[n_rc,~] = fn_num_subplots(numel(reg_lab));

for reg_ix = 1:numel(reg_lab)
    subplot(n_rc(1),n_rc(2),reg_ix); hold on;
    violin_data = struct;
    for cond_ix = 1:numel(cond_lab)
        violin_name = strrep(cond_names{cond_ix},' ','_');
        if any(cond_idx==cond_ix)
            violin_data.(violin_name) = model(cond_idx==cond_ix,reg_ix);
        else
            % Add empty (zero) data to prevent Violin error (e.g., neutral
            %   in early versions of the task
            violin_data.(violin_name) = 0;
        end
    end
    violins = violinplot(violin_data);
    
    % Plot properties
    legend_obj = cell(size(cond_lab));
    for cond_ix = 1:numel(cond_lab)
        % Change scatter colors to mark condition
        violins(cond_ix).ViolinColor = cond_colors{cond_ix};
        violins(cond_ix).BoxPlot.FaceColor = cond_colors{cond_ix};
        violins(cond_ix).EdgeColor = cond_colors{cond_ix};
        % Grab violin for legend
        legend_obj{cond_ix} = violins(cond_ix).ViolinPlot;
    end
    title(reg_names{reg_ix});
    % if strcmp(reg_lab{reg_ix},'EV'); legend([legend_obj{:}],cond_lab); end
    ax = gca;
    ax.XAxis.TickLabelInterpreter = 'none';
    ax.XAxis.TickLabelRotation = 45;
    ax.FontSize = 16;
end
saveas(gcf,[fig_dir fig_name '.png']);

%% Save Results
model_out_fname = [SBJ_vars.dirs.models SBJ '_model_' mdl.model_id '.mat'];
fprintf('Saving %s\n',model_out_fname);
if any(strcmp(reg_lab,'EV'))
    save(model_out_fname,'-v7.3','model','betas');
else
    save(model_out_fname,'-v7.3','model');
end

end
