function SBJ08a_HFA_actv(SBJ,pipeline_id,an_id,actv_win)
% Calculates high frequency activity, computes cluster-based statistics, and plots the results
% clear all; %close all;
% Set up paths
addpath('/home/knight/hoycw/PRJ_Error/scripts/');
addpath('/home/knight/hoycw/PRJ_Error/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults
if ischar(actv_win); actv_win = str2num(actv_win); end

%% Data Preparation
SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run /home/knight/hoycw/PRJ_Error/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
% eval(['run /home/knight/hoycw/PRJ_Error/scripts/proc_vars/' pipeline_id '_proc_vars.m']);

% Load Data
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',pipeline_id,'.mat'));
load(strcat(SBJ_vars.dirs.events,SBJ,'_trl_info_final.mat'));

%% Select Channel(s)
cfgs = [];
cfgs.channel = SBJ_vars.ch_lab.ROI;
roi = ft_selectdata(cfgs,data);

%% Cut into Trials
% Pad trial_lim_s by cfg_hfa.t_ftimwin/2 to avoid NaNs in epoch of interest
% Add 10 ms just because trimming back down to trial_lim_s exactly leave
% one NaN on the end, so smoothing will NaN out everything
trial_lim_s_pad = [trial_lim_s(1)-max(cfg_hfa.t_ftimwin)/2 trial_lim_s(2)+max(cfg_hfa.t_ftimwin)/2+0.01];

% Always normalize to pre-stimulus baseline for HFA
bsln_events = trl_info.trl_onset;
if strcmp(event_type,'stim')
    % Cut to desired trial_lim_s
    roi_trl = fn_ft_cut_trials_equal_len(roi,bsln_events,trl_info.cond_n,...
        round(trial_lim_s_pad*roi.fsample));
elseif strcmp(event_type,'resp')
    % Check that baseline will be included in trial_lim_s
    if trial_lim_s(1)>bsln_lim(1)
        error(['ERROR: trial_lim_s does not include bsln_lim for an_id = ' an_id]);
    end
    % Cut out to max_RT+trial_lim_s(2)+max(cfg_hfa.t_ftimwin)
    max_RT  = max(trl_info.rt);
    roi_trl = fn_ft_cut_trials_equal_len(roi,bsln_events,trl_info.cond_n,...
        round([trial_lim_s_pad(1) max_RT+trial_lim_s_pad(2)]*roi.fsample));
else
    error(['Unknown event_type: ' event_type]);
end

%% Compute HFA
fprintf('===================================================\n');
fprintf('------------------ HFA Calculations ---------------\n');
fprintf('===================================================\n');
if strcmp(HFA_type,'multiband')
    cfg_hfa.trials = 'all';
    hfa = ft_freqanalysis(cfg_hfa, roi_trl);
elseif strcmp(HFA_type,'broadband')
    error('Stop using broadband and use multitapers you dummy!');
    %         % Filter to HFA band
    %         cfgpp = [];
    %         cfgpp.hpfilter  = 'yes';
    %         cfgpp.hpfreq    = 70;
    %         cfgpp.lpfilter  = 'yes';
    %         cfgpp.lpfreq    = 150;
    %         roi = ft_preprocessing(cfgpp,roi);
    %         % Hilbert method to extract power
else
    error('Unknown HFA_type provided');
end

% Trim back down to original trial_lim_s to exclude NaNs
if strcmp(event_type,'stim')
    cfg_trim = [];
    cfg_trim.latency = trial_lim_s;
    hfa = ft_selectdata(cfg_trim,hfa);
elseif strcmp(event_type,'resp')
    cfg_trim = [];
    cfg_trim.latency = [trial_lim_s(1) max_RT+trial_lim_s(2)];
    hfa = ft_selectdata(cfg_trim,hfa);
else
    error(['Unknown event_type: ' event_type]);
end

%% Baseline Correction
fprintf('===================================================\n');
fprintf('---------------- Baseline Correction --------------\n');
fprintf('===================================================\n');
switch bsln_type
    case {'zboot', 'zscore'}
        hfa = fn_bsln_ft_tfr(hfa,bsln_lim,bsln_type,n_boots);
    case {'relchange', 'demean', 'my_relchange'}
        error(['bsln_type ' bsln_type ' is not compatible with one-sample t test bsln activation stats']);
%         cfgbsln = [];
%         cfgbsln.baseline     = bsln_lim;
%         cfgbsln.baselinetype = bsln_type;
%         cfgbsln.parameter    = 'powspctrm';
%         hfa = ft_freqbaseline(cfgbsln,hfa);
    otherwise
        error(['No baseline implemented for bsln_type: ' bsln_type]);
end

%% Smooth Power Time Series
if smooth_pow_ts
    fprintf('===================================================\n');
    fprintf('----------------- Filtering Power -----------------\n');
    fprintf('===================================================\n');
    for ch_ix = 1:numel(hfa.label)
        for f_ix = 1:numel(hfa.freq)
            if strcmp(lp_yn,'yes') && strcmp(hp_yn,'yes')
                hfa.powspctrm(:,ch_ix,f_ix,:) = fn_EEGlab_bandpass(...
                    hfa.powspctrm(:,ch_ix,f_ix,:), roi.fsample, hp_freq, lp_freq);
            elseif strcmp(lp_yn,'yes')
                hfa.powspctrm(:,ch_ix,f_ix,:) = fn_EEGlab_lowpass(...
                    squeeze(hfa.powspctrm(:,ch_ix,f_ix,:)), roi.fsample, lp_freq);
            elseif strcmp(hp_yn,'yes')
                error('Why are you only high passing?');
            else
                error('Why did you say yes smooth but no to both low and high pass?');
            end
        end
    end
end

%% Merge multiple bands
if strcmp(HFA_type,'multiband')
    cfg_avg = [];
    cfg_avg.freq = 'all';
    cfg_avg.avgoverfreq = 'yes';
    hfa = ft_selectdata(cfg_avg,hfa);
end

%% Re-align to event of interest if necessary (e.g., response)
if strcmp(event_type,'resp')
    hfa = fn_realign_tfr_s2r(hfa,trl_info.rt,stat_lim);
elseif ~strcmp(event_type,'stim')
    error(['ERROR: unknown event_type ' event_type]);
end

%% Run Statistics
fprintf('===================================================\n');
fprintf('--------------------- Statistics ------------------\n');
fprintf('===================================================\n');

% Select data in stat window
cfg_trim = [];
cfg_trim.latency = stat_lim;
hfa_stat = ft_selectdata(cfg_trim,hfa);

actv_ch = {};
actv_ch_epochs = {};
for ch_ix = 1:numel(hfa_stat.label)
    % Compute t-test per time point
    n_tests   = size(hfa_stat.powspctrm,4);
    pvals     = NaN([1 n_tests]);
    for time_ix = 1:n_tests
        [~, pvals(time_ix)] = ttest(squeeze(hfa_stat.powspctrm(:,ch_ix,1,time_ix)));
    end
    
    % Find epochs with significant task activations
%     [~, qvals] = mafdr(pvals); % Errors on some random ch (e.g., ROF8-9
%     in IR32), so I'm trying the below function
%     [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);
    [~, ~, ~, qvals] = fdr_bh(pvals);
    actv_mask = qvals<=0.05;
    actv_chunks = fn_find_chunks(actv_mask);
    actv_chunks(actv_mask(actv_chunks(:,1))==0,:) = [];
    actv_chunk_sz = diff(actv_chunks,1,2)+1;
    actv_epochs = actv_chunks(actv_chunk_sz > actv_win,:);
    if ~isempty(actv_epochs)
        actv_ch = {actv_ch{:} hfa_stat.label{ch_ix}};
        actv_ch_epochs = {actv_ch_epochs{:}, hfa_stat.time(actv_epochs)};
    end
end

%% Save Results
data_out_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_actv_ROI_',an_id,'_mn',num2str(actv_win),'.mat');
fprintf('===================================================\n');
fprintf('--- Saving %s ------------------\n',data_out_filename);
fprintf('===================================================\n');
save(data_out_filename,'-v7.3','hfa','actv_ch','actv_ch_epochs');

end
