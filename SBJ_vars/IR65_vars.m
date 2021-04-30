%% IR65 Processing Variables
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ        = 'IR65';
SBJ_vars.raw_file   = {'2017120714_0001.besa'};
SBJ_vars.block_name = {''};
SBJ_vars.low_srate  = [0];
SBJ_vars.log_fname  = {'865_response_log_20171207130759.txt'};

SBJ_vars.dirs.SBJ     = [root_dir 'PRJ_Error/data/' SBJ_vars.SBJ '/'];
SBJ_vars.dirs.raw     = [SBJ_vars.dirs.SBJ '00_raw/'];
SBJ_vars.dirs.import  = [SBJ_vars.dirs.SBJ '01_import/'];
SBJ_vars.dirs.preproc = [SBJ_vars.dirs.SBJ '02_preproc/'];
SBJ_vars.dirs.events  = [SBJ_vars.dirs.SBJ '03_events/'];
SBJ_vars.dirs.proc    = [SBJ_vars.dirs.SBJ '04_proc/'];
SBJ_vars.dirs.recon   = [SBJ_vars.dirs.SBJ '05_recon/'];
if ~exist(SBJ_vars.dirs.import,'dir')
    mkdir(SBJ_vars.dirs.import);
end
if ~exist(SBJ_vars.dirs.preproc,'dir')
    mkdir(SBJ_vars.dirs.preproc);
end
if ~exist(SBJ_vars.dirs.events,'dir')
    mkdir(SBJ_vars.dirs.events);
end
if ~exist(SBJ_vars.dirs.proc,'dir')
    mkdir(SBJ_vars.dirs.proc);
end
if ~exist(SBJ_vars.dirs.recon,'dir')
    mkdir(SBJ_vars.dirs.recon);
end

SBJ_vars.dirs.raw_filename = strcat(SBJ_vars.dirs.raw,SBJ_vars.raw_file);

%--------------------------------------
% Channel Selection
%--------------------------------------
hdr = ft_read_header(SBJ_vars.dirs.raw_filename);
SBJ_vars.orig_n_ch = length(hdr.label);
SBJ_vars.orig_n_samples = hdr.nSamples;
SBJ_vars.orig_srate = hdr.Fs;
clear hdr;

SBJ_vars.ch_lab.probes = {};
SBJ_vars.ch_lab.ref_type     = {};
SBJ_vars.ch_lab.ROI    = {};
SBJ_vars.ch_lab.eeg_ROI = {};

%SBJ_vars.ch_lab.prefix = ''; %'POL ';    % before every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.suffix = ''; %'-Ref';    % after every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.mislabel = {};% {{'ch1bad','ch1good'}};

SBJ_vars.ch_lab.ref_exclude = {}; %exclude from the CAR
SBJ_vars.ch_lab.bad = {...
    };
SBJ_vars.ch_lab.eeg = {};
% SBJ_vars.ch_lab.CZ_lap_ref = {};
SBJ_vars.ch_lab.eog = {};
SBJ_vars.ch_lab.photod = {};

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
SBJ_vars.analysis_time = {};
SBJ_vars.ignore_trials = [];
if numel(SBJ_vars.analysis_time) ~= numel(SBJ_vars.raw_file) || numel(SBJ_vars.raw_file) ~= numel(SBJ_vars.block_name)
    error('Mismatch number of runs to concatenate!');
end

%--------------------------------------
% Artifact Rejection Parameters
%--------------------------------------
% SBJ_vars.artifact_params.std_limit_raw = 7;
% SBJ_vars.artifact_params.hard_threshold_raw = 1000;

% SBJ_vars.artifact_params.std_limit_diff = 7;
% SBJ_vars.artifact_params.hard_threshold_diff = 100;

%--------------------------------------
% Trials to Reject
%--------------------------------------
% SBJ_vars.trial_reject_n = [];
