function [trl_info] = SBJ03_behav_parse(SBJ, pipeline_id, plot_it, save_it)
%
% SBJ [str] -- uniquely identifies the subject, e.g., 'IR54'
% pipeline_id [str] -- name of the processing pipeline (to get resample rate)
% plot_it [0/1] -- optional. plot_it = 1 to plot detected events
% save_it [0/1] -- whether to save output

min_event_length = 0.2; % photodiode elevation must be at least this long (in sec)
trig_dur_trl = 0.3;     % duration of trial onset trigger (in sec)
trig_dur_fdb = 0.8;     % duration of feedback onset trigger (in sec)
n_hist_bins = 100; % Set this to a higher number for more amplitudes. More bins can mean more errors though.

%% File paths
SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

utils_dir_name = '/home/knight/hoycw/PRJ_Error/scripts/utils/';
evnt_filename = [SBJ_vars.dirs.import SBJ '_evnt.mat'];
csv_filename = [SBJ_vars.dirs.events SBJ '_behav.csv'];
output_filename = [SBJ_vars.dirs.events SBJ '_trl_info_auto.mat'];

% Import helper functions to Matlab path
addpath(genpath(utils_dir_name));

%% Determine event onset sample points
% Load preprocessing variables to get sampling rate
proc_vars_cmd = ['run /home/knight/hoycw/PRJ_Error/scripts/proc_vars/' pipeline_id '_proc_vars.m'];
eval(proc_vars_cmd);
nrl_srate = proc_vars.resample_freq;

% Load input data
fprintf('Loading %s\n',evnt_filename);
load(evnt_filename);
n_samples = size(evnt.trial{1},2);
s_rate = evnt.fsample;
len_in_sec = size(evnt.trial{1},2)/evnt.fsample;
data_photo = evnt.trial{1};
data_photo_orig = data_photo;

% Bring data down to zero
data_photo = data_photo - min(data_photo);

% Read photodiode data
fprintf('\tReading photodiode data\n');
min_event_length = min_event_length * s_rate;    %trial must be at least 0.8 sec (actually ~1.5s?)
% Decimate to around 1000Hz (don't worry about aliasing)
if s_rate > 1000
  decimate_v = floor(s_rate/1000);
  data_photo_d = data_photo(1:decimate_v:end);
  min_event_length = floor(min_event_length/decimate_v);
else
    fprintf('================================================\n');
    fprintf('WARNING: evnt sample rate is 1000 or less!!\n');
    data_photo_d = data_photo;
    decimate_v = 1;
end

%% Process photodiode
[all_onsets, all_offsets, data_shades] = read_photodiode(data_photo_d, min_event_length, 2, plot_it);  %2 different shades (bsln, event), so set to 0
clear data_photo;

%% Separate out trial and feedback onsets
trl_onsets = zeros(size(all_onsets{2},1)/2,1);
fdb_onsets = zeros(size(all_onsets{2},1)/2,1);
trl_ix = 0;
fdb_ix = 0;
for evnt_ix = 1:numel(all_onsets{2})
    evnt_len = all_offsets{2}(evnt_ix)-all_onsets{2}(evnt_ix);
    if evnt_len < mean([trig_dur_trl trig_dur_fdb])*1000
        trl_ix = trl_ix+1;
%         fprintf('trl_ix %d, evnt_len = %d\n',trl_ix,evnt_len);
        trl_onsets(trl_ix) = all_onsets{2}(evnt_ix);
    else
        fdb_ix = fdb_ix+1;
%         fprintf('fdb_ix %d, evnt_len = %d\n',fdb_ix,evnt_len);
        fdb_onsets(fdb_ix) = all_onsets{2}(evnt_ix);
    end
end

% % Diff to get edges which correspond to onsets and offsets
% data_shades = [diff(data_shades) 0]; % Add a point because diff removes one
% trl_onsets = find(data_shades>0)'; % 1 to 2,3,4 is word onset. Transpose to make column vector
trl_onsets = trl_onsets*decimate_v; % convert back to ms
fdb_onsets = fdb_onsets*decimate_v; % convert back to ms
fprintf('\t\tFound %d events in photodiode channel\n', length(all_onsets));
fprintf('\t\tFound %d trials (%d feedbacks) in photodiode channel\n', sum(trl_onsets~=0), sum(fdb_onsets~=0));

% Plot intervals between trial and feedback onsets
if plot_it
    figure;
    histogram(fdb_onsets-trl_onsets,n_hist_bins);
    title([SBJ ': ' num2str(numel(trl_onsets)) ' feedback_delay intervals']);
end

%% Read in log file
% Open file, check total # trials
trl_info = fn_load_trl_info_csv(csv_filename,[]);
fprintf('\t\tFound %d trials in log file\n', numel(trl_info.trl_n));

% Reload, this time removing trials to ignore
fprintf('\t\tIgnoring %d trials\n', length(SBJ_vars.ignore_trials));
if ~isempty(SBJ_vars.ignore_trials)
    trl_info = fn_load_trl_info_csv(csv_filename,SBJ_vars.ignore_trials);
end

fprintf('\t\tKeeping %d trials from log file\n', numel(trl_info.trl_n));

%% Check if log and photodiode have different n_trials, plot and error out
if(length(trl_info.trl_n) ~= length(trl_onsets))
  % Plot photodiode data
  plot_photo = data_photo_orig - min(data_photo_orig);
  plot_photo = plot_photo / (max(plot_photo)-min(plot_photo));
  plot_photo = plot_photo + 0.25;
  plot(plot_photo, 'k'); hold on;
  % Plot word onsets
  for word_n = 1:length(trl_onsets)
    plot([trl_onsets(word_n) trl_onsets(word_n)],[1.30 1.40],'r','LineWidth',2);
    plot([trl_onsets(word_n) trl_onsets(word_n)],[-0.35 0.35],'r','LineWidth',2);
  end
  error('\nNumber of trials in log is different from number of trials found in event channel\n\n');
end

%% Put all the information into correct structures
fprintf('\t\tMean response time: %1.2f seconds from word onset\n', nanmean(trl_info.rt)); % Ignore NaNs
trl_info.trl_onset = trl_onsets;
trl_info.rsp_onset = trl_onsets + floor(trl_info.rt*s_rate);
trl_info.fdb_onset = fdb_onsets;

% Track no response trials
trl_info.rsp_onset(trl_info.rt<0) = -1;
trl_info.hit(trl_info.rt<0) = -1;
trl_info.score(trl_info.rt<0) = 0;

% Book keeping
trl_info.ignore_trials = SBJ_vars.ignore_trials;
trl_info.SBJ = SBJ;
trl_info.rt_type = 'auto_log';

%% Load Paradigm Parameters, Add to trl_info
prdm_vars = load([SBJ_vars.dirs.events SBJ '_prdm_vars.mat']);
prdm_fields = fieldnames(prdm_vars);
for f_ix = 1:numel(prdm_fields)
    trl_info.(prdm_fields{f_ix}) = prdm_vars.(prdm_fields{f_ix});
end

%% Convert condition codes from string to numbers
trl_info.cond_types = {'easy', 'hard'};
trl_info.cond_n = zeros(size(trl_info.trl_n));
for t_ix = 1:numel(trl_info.trl_n)
    trl_info.cond_n(t_ix) = find(strcmp(trl_info.cond_types,trl_info.cond{t_ix}));
end

%% Convert samples from event channel sample rate to ecog channel sample rate
s_rate_ratio = nrl_srate / evnt.fsample;
trl_info.trl_onset = round(trl_info.trl_onset*s_rate_ratio);
trl_info.rsp_onset = round(trl_info.rsp_onset*s_rate_ratio);
trl_info.fdb_onset = round(trl_info.fdb_onset*s_rate_ratio);
trl_info.sample_rate = nrl_srate;

%% Save results
if save_it
    save(output_filename, 'trl_info');
end

%% Plot results
if(plot_it ~= 0)
  figure;%('Position', [100 100 1200 800]);
  hold on;
  plot([0 len_in_sec],[0 0],'k');
  
  % Plot photodiode data
  plot_photo = data_photo_orig - min(data_photo_orig);
  plot_photo = plot_photo / (max(plot_photo)-min(plot_photo));
  plot_photo = plot_photo + 0.25;
  plot(linspace(0, len_in_sec, n_samples), plot_photo, 'Color', [0.5 0.8 0.8]);
  plot([0 len_in_sec],[0.25 0.25],'k');
  
  % Plot word onsets
  for word_n = 1:length(trl_info.trl_onset)
    plot([trl_info.trl_onset(word_n)/s_rate_ratio trl_info.trl_onset(word_n)/s_rate_ratio]/s_rate,[1.30 1.40],'b','LineWidth',2);
    plot([trl_info.trl_onset(word_n)/s_rate_ratio trl_info.trl_onset(word_n)/s_rate_ratio]/s_rate,[-0.35 0.35],'b','LineWidth',2);
  end
  
  % Plot resp onsets
  for resp_n = 1:length(trl_info.rsp_onset)
    plot([trl_info.rsp_onset(resp_n)/s_rate_ratio trl_info.rsp_onset(resp_n)/s_rate_ratio]/s_rate,[1.35 1.45],'g','LineWidth',2);
    plot([trl_info.rsp_onset(resp_n)/s_rate_ratio trl_info.rsp_onset(resp_n)/s_rate_ratio]/s_rate,[-0.30 0.30],'g','LineWidth',2);
  end
  
end


