function SBJ07a_HFA_save(SBJ,proc_id,an_id)
% Calculates high frequency activity, computes cluster-based statistics, and plots the results
% clear all; %close all;
% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Data Preparation
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);

% Load Data
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
load(strcat(SBJ_vars.dirs.events,SBJ,'_trl_info_final.mat'));

% Toss sampleinfo which can mess up trial cutting
if isfield(data,'sampleinfo')
    data = rmfield(data,'sampleinfo');
end

%% Select Channel(s)
cfgs = [];
cfgs.channel = SBJ_vars.ch_lab.ROI;
roi = ft_selectdata(cfgs,data);
clear data;

%% Cut into Trials
% Padding:
%   At a minimum, trial_lim_s must extend 1/2*max(filter_window) prior to 
%   the first data point to be estimated to avoid edges in the filtering window.
%   (ft_freqanalysis will return NaN for partially empty windows, e.g. an edge pre-trial,
%   but ft_preprocessing would return a filtered time series with an edge artifact.)
%   Also note this padding buffer should be at least 3x the slowest cycle
%   of interest.
if strcmp(HFA_type,'multiband')
    pad_len = 0.5*max(cfg_hfa.t_ftimwin)*3;
elseif any(strcmp(HFA_type,{'broadband','hilbert'}))
    % add 250 ms as a rule of thumb, or onger if necessary
    pad_len = 0.5*max([1/min(fois)*3 0.25]);
end
% Cut data to bsln_lim to be consistent across S and R locked (confirmed below)
%   Add extra 10 ms just because trimming back down to trial_lim_s exactly leave
%   one NaN on the end (smoothing that will NaN out everything)
trial_lim_s_pad = [min(bsln_lim)-pad_len trial_lim_s(2)+pad_len+0.01];

% Always normalize to pre-stimulus baseline for HFA
bsln_events = trl_info.trl_onset;
if strcmp(event_type,'stim')
    % Check that baseline will be included in data cut to trial_lim_s
    if trial_lim_s(1) < bsln_lim(1)
        error(['ERROR: trial_lim_s does not include bsln_lim for an_id = ' an_id]);
    end
    % Cut to desired trial_lim_s
    roi_trl = fn_ft_cut_trials_equal_len(roi,bsln_events,fn_condition_index('DifOut',trl_info),...
        round(trial_lim_s_pad*roi.fsample));
elseif strcmp(event_type,'resp')
    % Check that baseline will be included in data cut to trial_lim_s
    if trial_lim_s(1)+min(trial_info.response_time) < bsln_lim(1)
        error(['ERROR: trial_lim_s does not include bsln_lim for an_id = ' an_id]);
    end
    % Cut to max_RT+trial_lim_s(2) to include S baseline + full R-locked trial_lim_s
    max_RT  = max(trl_info.rt);
    roi_trl = fn_ft_cut_trials_equal_len(roi,bsln_events,fn_condition_index('DifOut',trl_info),...
        round([trial_lim_s_pad(1) max_RT+trial_lim_s_pad(2)]*roi.fsample));
elseif strcmp(event_type,'fb')
    % Cut from S baseline to full fb_locked:
    %   [trl_onset+trial_lim_s(1) to fb+trial_lim_s(2)]
    roi_trl = fn_ft_cut_trials_equal_len(roi,bsln_events,fn_condition_index('DifOut',trl_info),...
        round([trial_lim_s_pad(1) trl_info.prdm.trl_len+trial_lim_s_pad(2)]*roi.fsample));
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
elseif strcmp(HFA_type,'hilbert')
    % Create fake ft_freqanalysis struct
    hfa.label = roi_trl.label;
    hfa.freq  = fois;
    hfa.time  = roi_trl.time{1};
    hfa.powspctrm = zeros([numel(roi_trl.trial) numel(roi_trl.label) numel(fois) numel(roi_trl.time{1})]);
    hfa.dimord = 'rpt_chan_freq_time';
    hfa.trialinfo = roi_trl.trialinfo;
    for f_ix = 1:numel(fois)
        cfg_hfa.bpfreq = bp_lim(f_ix,:);
        cfg_hfa.hilbert = 'abs';
        fprintf('\n------> %s filtering: %.03f - %.03f\n', HFA_type, bp_lim(f_ix,1), bp_lim(f_ix,2));
        hfa_tmp = ft_preprocessing(cfg_hfa,roi_trl);
        for t_ix = 1:numel(roi_trl.trial)
            hfa.powspctrm(t_ix,:,f_ix,:) = hfa_tmp.trial{t_ix};
        end
    end
    clear hfa_tmp;
elseif strcmp(HFA_type,'broadband')
    error('Stop using broadband and account for 1/f you dummy!');
    %         % Filter to single HFA band
    %         cfgpp=[];.hpfilter='yes';.hpfreq=70;.lpfilter='yes';lpfreq=150;
    %         roi = ft_preprocessing(cfgpp,roi);
else
    error('Unknown HFA_type provided');
end

% Trim back down to original trial_lim_s to exclude NaNs
cfg_trim = [];
if strcmp(event_type,'stim')
    cfg_trim.latency = trial_lim_s;
elseif strcmp(event_type,'resp') && strcmp(bsln_evnt,'stim')
    cfg_trim.latency = [bsln_lim(1) max_RT+trial_lim_s(2)];
elseif strcmp(event_type,'fb') && strcmp(bsln_evnt,'stim')
    cfg_trim.latency = [bsln_lim(1) trl_info.prdm.trl_len+trial_lim_s(2)];
else
    error('mismatched event without S-locked baseline!');
end
hfa = ft_selectdata(cfg_trim,hfa);

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
    % error catches
    if ~strcmp(lp_yn,'yes')
        if strcmp(hp_yn,'yes')
            error('Why are you only high passing?');
        else
            error('Why is smooth_pow_ts yes but no lp or hp?');
        end
    end
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
                    hfa.powspctrm(:,ch_ix,f_ix,:), roi.fsample, lp_freq);
            else
                error('weird non-Y/N filtering options!');
            end
        end
    end
end

%% Merge multiple bands
cfg_avg = [];
cfg_avg.freq = 'all';
cfg_avg.avgoverfreq = 'yes';
hfa = ft_selectdata(cfg_avg,hfa);

%% Re-align to event of interest if necessary (e.g., response)
if strcmp(event_type,'resp')
    hfa = fn_realign_tfr_s2r(hfa,trl_info.rt,trial_lim_s);
elseif strcmp(event_type,'fb')
    hfa = fn_realign_tfr_s2r(hfa,...
        ones(size(trl_info.rt))*trl_info.prdm.target+trl_info.prdm.fb_delay,...
        trial_lim_s);
elseif ~strcmp(event_type,'stim')
    error(['ERROR: unknown event_type ' event_type]);
end

%% Downsample
if resample_ts && hfa.fsample~=resample_freq
    cfgrs = [];
    cfgrs.resamplefs = resample_freq;
    cfgrs.detrend = 'no';
    hfa = ft_resampledata(cfgrs, hfa);
end

%% Save Results
data_out_fname = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'.mat');
fprintf('===================================================\n');
fprintf('--- Saving %s ------------------\n',data_out_fname);
fprintf('===================================================\n');
save(data_out_fname,'-v7.3','hfa');

end
