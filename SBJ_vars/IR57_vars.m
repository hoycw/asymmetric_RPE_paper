%% ${SBJ} Processing Variables
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ        = 'IR57';
SBJ_vars.raw_file   = {'2017032208_0002.besa'};
SBJ_vars.block_name = {''};
SBJ_vars.low_srate  = [0];
SBJ_vars.log_fname  = {'857_response_log_20170322112243_CWHedit.txt'};

SBJ_vars.dirs.SBJ     = ['/home/knight/hoycw/PRJ_Error/data/' SBJ_vars.SBJ '/'];
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

SBJ_vars.ch_lab.probes = {'RSM','RAC','ROF','RIN','RTI','RAM','RHH','RTH',...
                         'LSMA','LAC','LOF','LIN','LTI','LAM','LTH'};%'LHH' doesn't count because all elecs are bad
SBJ_vars.ref_types     = {'BP','BP','BP','BP','BP','BP','BP','BP',...
                          'BP','BP','BP','BP','BP','BP','BP'};
SBJ_vars.ch_lab.ROI    = {'RSM*','RAC*','ROF*','RIN*','RTI*',...
                              'LAM4-5','LAM5-6','LAC*','LOF*'};%LAM4,5 are inferior anterior insula
% PRJ_Stroop varaince rejected: {'-LOF1-2','-RTI2-3','-RIN4-5'};

%SBJ_vars.ch_lab.prefix = ''; %'POL ';    % before every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.suffix = ''; %'-Ref';    % after every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.mislabel = {};

% Based on below, LSMA1, RSM1, and RHH5 are left out (no pairs)
SBJ_vars.ch_lab.ref_exclude = {}; %exclude from the CAR
SBJ_vars.ch_lab.bad = {...
    'RHH1','RHH2','RHH3','RHH4','ROF1','ROF2','RAM1','RTH1','RTH2',...%epileptic
    'LHH1','LHH2','LHH3','LHH4','LHH5','LHH6','LHH7','LHH8','LHH9','LHH10','LOF1',...%epileptic
    'LTH1','LTH2','LTH3','LAM1','LAM2','LAM3',...%epileptic
    'RHH6','RHH7',...%bad line noise, RSM2-3 (real bad) and LSMA3 PSDs looks worse than surrounding, but maybe not so bad it can't be saved
    'RSM9','RSM10','RAC10','RAM10','RTH9','RTH10','LOF10','LAM10',...%out of brain
    'DC01','DC03','DC04','XREF','E',...% not real data
    'EKG'...
    };
SBJ_vars.ch_lab.eeg = {'FPZ','CZ','OZ','C3','C4','Z','FP1','FP2','T3','T4','O1','O2'};
SBJ_vars.ch_lab.eeg_bad = {'Z','T4','O2','T3'}; %take these out of the above .eeg field!
SBJ_vars.ch_lab.eog = {'LUC','LLC','RUC','RLC'};
SBJ_vars.ch_lab.photod = {'DC02'};

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
% most have only regular harmonics with normal width, and 120 and 240 are weak
% RBT, RPIN, RSMA,LAC have an extra peak at 200
% RHH6 has really bad at all harmonics, like LUE and other nonsense chan
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
% data starts ~62s, goes to ~1410
SBJ_vars.analysis_time = {[52 1420]};
SBJ_vars.ignore_trials = [];

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
% These should be indices AFTER SBJ05 has run!
% SBJ_vars.trial_reject_ix = [];
