%% SBJ00_cleaning_prep
% Load raw data, extract channel labels, downsample, bandstop line noise,
% and save the result for visual cleaning inspection.
function SBJ00_cleaning_prep(SBJ,plot_psd)
% Inputs:
%   SBJ [str]- name of the subject
%   plot_psd [str]- name of function to plot+save PSDs of raw channels
%       'all'- plots all PSDs on one figure
%       '1by1'- plots each PSD on a separate figure
%       any other string- nothing is plotted or saved
% Outputs:
%   SBJ_preclean.mat- the resampled and bandstop filtered data for visual cleaning
%   SBJ_raw_labels.mat- the raw labels for all channels, to inform channel selection
%   raw_psds/- PSDs of channels as determined by plot_psd

addpath('/home/knight/hoycw/PRJ_Error/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

% plot_psd      = '1by1';         % 'all','1by1','none'
psd_fig_type  = 'jpg';
resamp_it     = 1;
resample_freq = 1000;
filter_it     = 1;
notch_freqs   = [60 120 180 240 300];
bs_width      = 2;
bs_freq_lim   = NaN([length(notch_freqs) 2]);
for f_ix = 1:length(notch_freqs)
    bs_freq_lim(f_ix,:) = fn_freq_lim_from_CFBW(notch_freqs(f_ix),bs_width);
end

%% Processing
SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
fprintf('============== Processing %s, %s ==============\n',SBJ,SBJ_vars.raw_file);

%% Set up directories
% Assume only one block for now...
%         if length(SBJ_files{SBJ_ix}{2})>1
%             block_prefix = strcat('_B',num2str(file_ix));
psd_dir = strcat(SBJ_vars.dirs.import,'raw_psds/');
if ~exist(psd_dir,'dir')
    mkdir(psd_dir);
end

%% Load the data
if strcmp(SBJ_vars.raw_file(end-2:end),'mat')
    load(SBJ_vars.dirs.raw_filename);
else
    cfg = [];
    cfg.dataset = SBJ_vars.dirs.raw_filename;
    cfg.continuous = 'yes';
    cfg.channel = 'all';
    data = ft_preprocessing(cfg);
end

%% Check for problematic labels
for ch_ix = 1:numel(data.label)
    if ~isempty(find(data.label{ch_ix}=='/',1))   % These will mess up paths when saved
        data.label{ch_ix}(find(data.label{ch_ix}=='/')) = '_';
    end
end

%% Check noise profile
if strcmp(plot_psd,'all')
    fprintf('============== Plotting PSDs %s, %s ==============\n',SBJ,SBJ_vars.raw_file);
    fn_plot_PSD_all_save(data.trial{1},data.label,data.fsample,...
        strcat(psd_dir,SBJ,'_raw_psd',SBJ_vars.block_prefix,'_all.',psd_fig_type));
elseif strcmp(plot_psd,'1by1')
    fprintf('============== Plotting PSDs %s, %s ==============\n',SBJ,SBJ_vars.raw_file);
    fn_plot_PSD_1by1_save(data.trial{1},data.label,data.fsample,...
        strcat(psd_dir,SBJ,'_raw_psd',SBJ_vars.block_prefix),psd_fig_type);
end

%% Resample data to speed things up
if (resamp_it) && (data.fsample ~= resample_freq)
    cfg_resamp = [];
    cfg_resamp.resamplefs = resample_freq;
    cfg_resamp.detrend = 'yes';
    data = ft_resampledata(cfg_resamp, data);
end

%% Filter data for ease of viewing
fprintf('============== Bandstop filtering %s, %s ==============\n',SBJ,SBJ_vars.raw_file);
if (filter_it)% && (data_resamp.cfg.dftfreq ~= notch_freqs)
    cfg_bs = [];
    cfg_bs.continuous = 'yes';
    cfg_bs.bsfilter = 'yes';
    cfg_bs.bsfreq = bs_freq_lim;
    cfg_bs.demean = 'yes';
    data = ft_preprocessing(cfg_bs, data);
end

%% Save data
% Save preclean data
out_filename = strcat(SBJ_vars.dirs.preproc,SBJ,'_preclean',SBJ_vars.block_prefix,'.mat');
fprintf('============== Saving %s, %s ==============\n',out_filename,SBJ_vars.raw_file);
save(out_filename, '-v7.3', 'data');

% Save data labels
raw_labels = data.label;
label_filename = strcat(SBJ_vars.dirs.import,SBJ,'_raw_labels',SBJ_vars.block_prefix,'.mat');
% if raw_labels exit, load, compare, save if different
if ~exist(label_filename)
    save(label_filename,'raw_labels');
else
    error(strcat('ERROR: raw_labels filename already exists: ',label_filename));
end

end
