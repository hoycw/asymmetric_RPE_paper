%% CP24 Processing Variables
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
if isempty(strfind(path,'fieldtrip'))
    addpath(ft_dir);
    ft_defaults
end

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ        = 'CP24';
SBJ_vars.raw_file   = {'CP24_Dec6_error_R1.mat','CP24_Dec9_error_R2.mat'};
SBJ_vars.block_name = {'R1','R2'};
SBJ_vars.low_srate  = [0,0];
SBJ_vars.log_file   = {'224_response_log_20171206121023.txt','cp24_2_response_log_20171209120902.txt'};

SBJ_vars.dirs.SBJ     = [root_dir 'PRJ_Stroop/data/' SBJ_vars.SBJ '/'];
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

SBJ_vars.recon.surf_l     = [SBJ_vars.dirs.recon 'Surfaces/' SBJ_vars.SBJ '_cortex_lh.mat'];
SBJ_vars.recon.surf_r     = [SBJ_vars.dirs.recon 'Surfaces/' SBJ_vars.SBJ '_cortex_rh.mat'];
SBJ_vars.recon.elec_pat   = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_acpc_r.mat'];
SBJ_vars.recon.elec_mni_v = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_mni_v.mat'];
SBJ_vars.recon.elec_mni_s = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_mni_s.mat'];
SBJ_vars.recon.fs_T1      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_T1.mgz'];
SBJ_vars.recon.fs_DK      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_aparc+aseg.mgz'];
SBJ_vars.recon.fs_Dx      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_aparc.a2009s+aseg.mgz'];

%--------------------------------------
% Channel Selection
%--------------------------------------
SBJ_vars.ch_lab.probes     = {'RMT','RTO','RIHA','RIHP','ROF','RLF','LMT','LTO','LIHA','LIHP','LOF'};%,'LLFP','LLF'};
SBJ_vars.ch_lab.probe_type = {'ecog','ecog','ecog','ecog','ecog','ecog','ecog','ecog','ecog','ecog','ecog'};
SBJ_vars.ch_lab.ref_type   = {'CARall','CARall','CARall','CARall','CARall','CARall',...
                                'CARall','CARall','CARall','CARall','CARall'};
%SBJ_vars.ch_lab.ref_type   = {'CAR','CAR','CAR','CAR','CAR','CAR','CAR','CAR','CAR','CAR','CAR'};%,'CAR','CAR'};
if ~all(numel(SBJ_vars.ch_lab.probes)==[numel(SBJ_vars.ch_lab.probe_type) numel(SBJ_vars.ch_lab.ref_type)]); error('probes ~= type+ref');end;
SBJ_vars.ch_lab.ROI        = {'all'};%'RIHA*','RIHP*','ROF*','RLF*','LIHA*','LIHP*','LOF*'};
SBJ_vars.ch_lab.eeg_ROI    = {};

%SBJ_vars.ch_lab.prefix = 'POL ';    % before every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.suffix = '-Ref';    % after every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.mislabel = {{'RLT12','FPG12'},{'IH;L8','IHL8'}};

SBJ_vars.ch_lab.ref_exclude = {... %exclude from the CAR
     'RLF1',...% a litle spread from RLF2 sometimes
     'RTO3',...% still correlations about like EOG-ish channels
     'ROF1','ROF2','ROF3','LOF2','LOF3','LOF4',...% EOG big stuff, weak correlations
     'RTO6','RTO7',...% HF noise
     'RTO8'...% occassional spikes
     }; %exclude from the CAR
% emodim ref_exclude:
%     'RTO3',...% still correlations about like EOG-ish channels
%     'ROF1','ROF2','ROF3','LOF2','LOF3','LOF4',...% EOG big stuff, weak correlations
%     'RTO6','RTO7',...% HF noise
%     'RTO8'...% occassional spikes
%     }; %exclude from the CAR
SBJ_vars.ch_lab.bad = {...
    'LTO2',...% spiking
    'RLF2','RTO4',...% epileptic
    'RTO1','RTO2',...% flat channels, no signal
    'LLF*','LLFP*',...% noisy (on second amplifier)
    'DC01','DC04','E','EEG Mark1','EEG Mark2','-','Events/Markers',...%not real data
    'Events_Markers','Mark1','Mark2','X54','Events'...% Not real data
    };
% only bad channel in original bad: 'RLF2'...% epileptic
% LLFP and LLF have massive noise, maybe can't save them
% emodim .bad:
%    'LTO2',...% spiking
%    'RLF1','RLF2','RTO4',...% epileptic
%    'LLF*','LLFP*',...% noisy (on second amplifier)
%    'RTO1','RTO2',...% flat channels, no signal
% bad_codes: 1 = toss (epileptic or bad); 2 = suspicious; 3 = out of brain; 0 = junk
SBJ_vars.ch_lab.bad_type = {'bad','sus','out'};
SBJ_vars.ch_lab.bad_code = [1 1 1 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
if numel(SBJ_vars.ch_lab.bad)~=numel(SBJ_vars.ch_lab.bad_code);error('bad ~= bad_code');end
SBJ_vars.ch_lab.eeg = {};
% SBJ_vars.ch_lab.CZ_lap_ref = {};
SBJ_vars.ch_lab.eog = {};
SBJ_vars.ch_lab.photod = {'DC02'};

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
% see cut scripts for timing info
SBJ_vars.analysis_time = {{[1 1289]},{[1 1429]}};

%--------------------------------------
% Artifact Rejection Parameters
%--------------------------------------
SBJ_vars.artifact_params.std_limit_raw = 7;
SBJ_vars.artifact_params.hard_threshold_raw = 1000;

SBJ_vars.artifact_params.std_limit_diff = 7;
SBJ_vars.artifact_params.hard_threshold_diff = 100;

%--------------------------------------
% Trials to Reject
%--------------------------------------
SBJ_vars.trial_reject_n = [];
