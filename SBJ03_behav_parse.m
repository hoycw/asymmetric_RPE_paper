function [trl_info] = SBJ03_behav_parse(SBJ, block, proc_id, plot_it, save_it)
%
% SBJ [str] -- uniquely identifies the subject, e.g., 'IR54'
% block [int] -- index of which block of data should be analyzed
% proc_id [str] -- name of the processing pipeline (to get resample rate)
% plot_it [0/1] -- optional. plot_it = 1 to plot detected events
% save_it [0/1] -- whether to save output

n_hist_bins = 100; % Set this to a higher number for more amplitudes. More bins can mean more errors though.

%% File paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
if numel(SBJ_vars.raw_file)>1
    block_suffix = strcat('_',SBJ_vars.block_name{block});
else
    block_suffix = SBJ_vars.block_name{block};   % should just be ''
end

evnt_fname   = [SBJ_vars.dirs.preproc SBJ '_evnt_clean',block_suffix,'.mat'];
bhv_fname    = [SBJ_vars.dirs.events SBJ '_behav',block_suffix,'.csv'];
output_fname = [SBJ_vars.dirs.events SBJ '_trl_info_auto',block_suffix,'.mat'];

% Import helper functions to Matlab path
addpath(genpath([root_dir 'PRJ_Error/scripts/utils/']));

%% Determine event onset sample points
% Load preprocessing variables to get sampling rate
eval(['run ' root_dir 'PRJ_Error/scripts/proc_vars/' proc_id '_proc_vars.m']);
if any(SBJ_vars.low_srate)
    nrl_srate = SBJ_vars.low_srate(block);
else
    nrl_srate = proc_vars.resample_freq;
end

% Load paradigm variables
prdm_vars = load([SBJ_vars.dirs.events SBJ '_prdm_vars.mat']);
if isfield(prdm_vars,'trig_dur')
    trig_dur_trl = prdm_vars.trig_dur;
else
    trig_dur_trl = 0.2;
end
trig_dur_fb = prdm_vars.fb;

% Load input data
fprintf('Loading %s\n',evnt_fname);
load(evnt_fname);
n_samples  = size(evnt.trial{1},2);
evnt_srate = evnt.fsample;
len_in_sec = size(evnt.trial{1},2)/evnt_srate;
data_photo = evnt.trial{1};
data_photo_orig = data_photo;

% Bring data down to zero
data_photo = data_photo - min(data_photo);

% Read photodiode data
fprintf('\tReading photodiode data\n');
% photodiode elevation must be at least this long (in sec)
%   assumes trial onset trigger = 0.2
min_event_length = 0.1*evnt_srate;
if evnt_srate < 1000
    fprintf('================================================\n');
    fprintf('WARNING: evnt sample rate is less than 1000 Hz!! s_rate = %f\n',evnt_srate);
end

%% Process photodiode
% 2 different shades (bsln, event), so set to 0
[all_onsets, all_offsets, ~] = read_photodiode(data_photo, min_event_length, 2, plot_it);
if save_it
    fig_fname = [SBJ_vars.dirs.events SBJ '_photo_segmentation' block_suffix '.fig'];
    saveas(gcf,fig_fname);
end
clear data_photo;

%% Separate out trial and feedback onsets
trl_onsets = zeros(size(all_onsets{2},1)/2,1);
fb_onsets  = zeros(size(all_onsets{2},1)/2,1);
trl_ix = 0;
fb_ix  = 0;
for evnt_ix = 1:numel(all_onsets{2})
    evnt_len = all_offsets{2}(evnt_ix)-all_onsets{2}(evnt_ix);
    if evnt_len < mean([trig_dur_trl trig_dur_fb])*1000
        trl_ix = trl_ix+1;
%         fprintf('trl_ix %d, evnt_len = %d\n',trl_ix,evnt_len);
        trl_onsets(trl_ix) = all_onsets{2}(evnt_ix);
    else
        fb_ix = fb_ix+1;
%         fprintf('fb_ix %d, evnt_len = %d\n',fb_ix,evnt_len);
        fb_onsets(fb_ix) = all_onsets{2}(evnt_ix);
    end
end

% % Diff to get edges which correspond to onsets and offsets
% data_shades = [diff(data_shades) 0]; % Add a point because diff removes one
% trl_onsets = find(data_shades>0)'; % 1 to 2,3,4 is word onset. Transpose to make column vector
fprintf('\t\tFound %d events in photodiode channel\n', length(all_onsets));
fprintf('\t\tFound %d trials (%d feedbacks) in photodiode channel\n', sum(trl_onsets~=0), sum(fb_onsets~=0));

% Plot intervals between trial and feedback onsets
if plot_it
    figure;
    histogram((fb_onsets-trl_onsets)/evnt_srate,n_hist_bins);
    title([SBJ ': ' num2str(numel(trl_onsets)) ' feedback delay intervals']);
end

%% Read in log file
trl_info = fn_load_behav_csv(bhv_fname,ignore_trials);
fprintf('\t\tKeeping %d trials from log file\n', numel(trl_info.trl_n));

%% Check if log and photodiode have different n_trials, plot and error out
if(length(trl_info.trl_n) ~= length(trl_onsets))
    figure;
    % Plot photodiode data
    plot_photo = data_photo_orig - min(data_photo_orig);
    plot_photo = plot_photo / (max(plot_photo)-min(plot_photo));
    plot_photo = plot_photo + 0.25;
    plot(plot_photo, 'k'); hold on;
    % Plot word onsets
    for trial_ix = 1:length(trl_onsets)
        plot([trl_onsets(trial_ix) trl_onsets(trial_ix)],[1.30 1.40],'r','LineWidth',2);
        plot([trl_onsets(trial_ix) trl_onsets(trial_ix)],[-0.35 0.35],'r','LineWidth',2);
    end
    error('\nNumber of trials in log is different from number of trials found in event channel\n\n');
end

%% Put all the information into correct structures
fprintf('\t\tMean response time: %1.2f seconds from word onset\n', nanmean(trl_info.rt)); % Ignore NaNs
trl_info.trl_onset = trl_onsets;
trl_info.rsp_onset = trl_onsets + floor(trl_info.rt*evnt_srate);
trl_info.fb_onset  = fb_onsets;

% Track no response trials
trl_info.rsp_onset(trl_info.rt<0) = -1;
trl_info.hit(trl_info.rt<0)       = -1;
trl_info.score(trl_info.rt<0)     = 0;

% Book keeping
trl_info.ignore_trials = ignore_trials;
trl_info.SBJ     = SBJ;
trl_info.rt_type = 'auto_log';
trl_info.prdm    = prdm_vars;

%% Convert samples from event channel sample rate to iEEG channel sample rate
s_rate_ratio = nrl_srate / evnt.fsample;
trl_info.trl_onset = round(trl_info.trl_onset*s_rate_ratio);
trl_info.rsp_onset = round(trl_info.rsp_onset*s_rate_ratio);
trl_info.fb_onset  = round(trl_info.fb_onset*s_rate_ratio);
trl_info.sample_rate = nrl_srate;

%% Save results
if save_it
    save(output_fname, '-v7.3', 'trl_info');
end

%% Plot results
if (plot_it ~= 0)
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
    for trial_ix = 1:length(trl_info.trl_onset)
        plot([trl_info.trl_onset(trial_ix)/s_rate_ratio trl_info.trl_onset(trial_ix)/s_rate_ratio]/evnt_srate,[1.30 1.40],'b','LineWidth',2);
        plot([trl_info.trl_onset(trial_ix)/s_rate_ratio trl_info.trl_onset(trial_ix)/s_rate_ratio]/evnt_srate,[-0.35 0.35],'b','LineWidth',2);
    end
    
    % Plot resp onsets
    for resp_n = 1:length(trl_info.rsp_onset)
        plot([trl_info.rsp_onset(resp_n)/s_rate_ratio trl_info.rsp_onset(resp_n)/s_rate_ratio]/evnt_srate,[1.35 1.45],'g','LineWidth',2);
        plot([trl_info.rsp_onset(resp_n)/s_rate_ratio trl_info.rsp_onset(resp_n)/s_rate_ratio]/evnt_srate,[-0.30 0.30],'g','LineWidth',2);
    end
    
    % Plot feedback onsets
    for fb_n = 1:length(trl_info.fb_onset)
        plot([trl_info.fb_onset(fb_n)/s_rate_ratio trl_info.fb_onset(fb_n)/s_rate_ratio]/evnt_srate,[1.4 1.5],'r','LineWidth',2);
        plot([trl_info.fb_onset(fb_n)/s_rate_ratio trl_info.fb_onset(fb_n)/s_rate_ratio]/evnt_srate,[-0.30 0.30],'r','LineWidth',2);
    end
    
    if save_it
        fig_fname = [SBJ_vars.dirs.events SBJ '_events' block_suffix '.fig'];
        saveas(gcf,fig_fname);
    end
end

end


