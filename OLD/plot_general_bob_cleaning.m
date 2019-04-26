%% Plot preprocessed data for cleaning with Bob
%   1) Load preproc data (downsample, highpass, trim, notch)
%   2) Plot data to identify bad channels and epochs
clear all; %close all;

%addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/_TOOLBOXES/'));
addpath('/home/knight/hoycw/Apps/fieldtrip/fieldtrip/');
ft_defaults

%%
SBJ   = 'IR48';
task  = 'Color_Wheel';
block = 1;
SBJ_dir = fullfile('/home/knight/hoycw/Data_Extraction/data/',SBJ,task);
preproc_dir = fullfile(SBJ_dir,'02_preproc/');
preclean_filename = strcat(preproc_dir,SBJ,'_',task,'_preclean_B',num2str(block),'.mat');

%% Load the data
fprintf('============== Loading %s, %s, %i ==============\n',SBJ, task, block);
load(preclean_filename);
data_notch = [];
data_notch.label = header_ecog.channel_labels;
data_notch.fsample = header_ecog.sample_rate;
data_notch.trial{1} = data_ecog;
data_notch.time{1} = linspace(0, size(data_ecog,2), size(data_ecog,2)); % ms time spacing

% cfg_notch = [];
% cfg_notch.continuous = 'yes';
% cfg_notch.dftfilter = 'yes'; % line noise removal using discrete fourier transform
% cfg_notch.dftfreq = [100 200 240];
% cfg_notch.demean = 'yes';
% data_notch = ft_preprocessing(cfg_notch, data_notch);

%% Plot data: check good channels, refs, artifacts, etc.
cfg_plot = [];
cfg_plot.viewmode   = 'vertical'; %trace per line
cfg_plot.continuous = 'yes';    % not trials
cfg_plot.plotlabels = 'yes';    %to know what to cut
cfg_plot.blocksize  = 10;        % good for bob, can adjust
cfg_plot.ylim       = [-100 100];

plot_ch = data_notch.label;

plot_out = ft_databrowser(cfg_plot, data_notch); %pulls up your plotted data

%% Set up files to prepare for KLA scripts
analysis_channels = 1:length(data_notch.label);
if strcmp(raw_filename(end-4:end),'.besa')
    data_ecog = data_notch.trial{:}(analysis_channels,:); %only analysis_channels
    %   BUT! analysis_channels variable should be original chan #s
    header_ecog.sample_rate = data_notch.fsample;
    header_ecog.orig_n_channels = size(data_notch.label,1);
    header_ecog.n_samples = size(data_ecog,2);
    header_ecog.length_in_seconds = size(data_notch.trial{:},2)/data_notch.fsample;
    header_ecog.original_sample_rate = data_notch.fsample;
    header_ecog.channel_labels = {data_notch.label{analysis_channels}};
elseif strcmp(raw_filename(end-3:end),'.edf')
    data_ecog = data_ecog(analysis_channels,:);
    header_ecog.channel_labels = header_ecog.channel_labels(analysis_channels);
end
header_ecog.n_channels = size(analysis_channels,2);
header_ecog.orig_channel_n = analysis_channels;
header_ecog.raw_file = raw_filename;
header_ecog.line_noise_freqs = notch_freqs;
if resample_it == 1
    header_ecog.original_sample_rate = orig_srate;
    header_ecog.sample_rate = resample_freq;
end

%% Save out those files
out_filename = strcat(out_dir,SBJ,'_',task,'_bob_clean_B',num2str(block),'.mat');
fprintf('============== Saving %s, %s, %i ==============\n',out_filename);
save(out_filename, 'plot_out', 'header_ecog', 'data_ecog');%, 'header_evnt', 'data_evnt');

clear data_ecog header_ecog data_raw data_notch data_resamp

%% OLD STUFF
% % Grab microphone recording separately to export and listen
% mic_data = data_evnt(mic_channel_n,:);
% %rescale to prevent clipping, add 0.05 fudge factor
% mic_data_rescale = mic_data./(max(abs(mic_data))+0.05);

% if strcmp(raw_filename(end-4:end),'.besa')
%     data_evnt = data_raw.trial{:}(event_channels,:);
%     header_evnt.sample_rate = data_raw.fsample;
%     header_evnt.orig_n_channels = size(data_raw.label,1);
%     header_evnt.n_samples = size(data_ecog,2);
%     header_evnt.length_in_seconds = size(data_raw.trial{:},2)/data_raw.fsample;
%     header_evnt.original_sample_rate = data_raw.fsample;
%     header_evnt.channel_labels = {data_raw.label{event_channels}};
% elseif strcmp(raw_filename(end-3:end),'.edf')
%     data_evnt = data_ecog(event_channels,:);
%     header_evnt = header_ecog;
%     header_evnt.channel_labels = header_ecog.channel_labels(event_channels);
% end
% header_evnt.n_channels = size(event_channels,2);
% header_evnt.orig_channel_n = event_channels;
% header_evnt.raw_file = raw_filename;

% data_out_filename = strcat(out_dir,raw_data_id,'.mat');
% proc_out_filename = strcat(out_dir,raw_data_id,'_proc_vars.mat');
% mic_data_filename = strcat(out_dir,raw_data_id,'_mic_recording.wav');

% save(proc_out_filename, 'analysis_time', 'analysis_channels', ...
%     'event_channels', 'reref_channels', 'reref_weights', 'photo_channel_n', 'mic_channel_n');
% audiowrite(mic_data_filename,mic_data_rescale,header_evnt.sample_rate);



% Call to extract_nk_...:
% subject_id, input_file_name, save_file_name, analysis_channels, event_channels, ...
%   'data_dir_name', ['../' subject_id], ...
%   'resample_rate', 1000, ...
%   'reref_channels', reref_channels, ...
%   'reref_name', 'WM', ...
%   'reref_weights', reref_weights, ...
%   'line_noise_freqs', line_noise_freqs, ...
%   'analysis_time', analysis_time, ...
%   'highpass_freq', hp_cutoff ...

% Contents of header that extract_nk_... produces:
%    header.recording_startdate  - [string], DD.MM.YY of original recording
%    header.recording_starttime  - [string], HH.MM.SS of original recording
%    header.orig_n_channels      - [integer], The number of channels originally recorded
%    header.n_channels           - [integer], The number of channels currently extracted
%    header.n_samples            - [integer], The number of samples per channel
%    header.length_in_seconds    - [float], The length of the extracted data in seconds
%    header.sample_rate          - [integer], The current sampling rate in Hz.
%    header.original_sample_rate - [integer], The recording sampling rate in Hz. Usually 5000.
%    header.channel_labels       - [cell array of strings], The labels of each extracted channel.
%    header.orig_channel_n       - [array of integers], The corresponding channel numbers of the
%                                  extracted data in the original data.

