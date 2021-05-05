%% IR74 Processing Variables
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
if isempty(strfind(path,'fieldtrip')); addpath(ft_dir); ft_defaults; end

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ        = 'IR74';
SBJ_vars.raw_file   = {'2008123123_0001.besa'};
SBJ_vars.block_name = {''};
SBJ_vars.low_srate  = [0];
SBJ_vars.log_fname  = {'ir742_response_log_20180327170333.txt'};

SBJ_vars.dirs.SBJ     = [root_dir 'PRJ_Error/data/' SBJ_vars.SBJ '/'];
SBJ_vars.dirs.raw     = [SBJ_vars.dirs.SBJ '00_raw/'];
SBJ_vars.dirs.import  = [SBJ_vars.dirs.SBJ '01_import/'];
SBJ_vars.dirs.preproc = [SBJ_vars.dirs.SBJ '02_preproc/'];
SBJ_vars.dirs.events  = [SBJ_vars.dirs.SBJ '03_events/'];
SBJ_vars.dirs.proc    = [SBJ_vars.dirs.SBJ '04_proc/'];
SBJ_vars.dirs.recon   = [SBJ_vars.dirs.SBJ '05_recon/'];
SBJ_vars.dirs.models  = [SBJ_vars.dirs.SBJ '06_models/'];
SBJ_vars.dirs.stats   = [SBJ_vars.dirs.SBJ '07_stats/'];
dirs_fields = fieldnames(SBJ_vars.dirs);
for field_ix = 1:numel(dirs_fields)
    if ~strcmp(dirs_fields{field_ix},'raw_filename') && ~exist(SBJ_vars.dirs.(dirs_fields{field_ix}),'dir')
        mkdir(SBJ_vars.dirs.(dirs_fields{field_ix}));
    end
end

SBJ_vars.dirs.raw_filename = strcat(SBJ_vars.dirs.raw,SBJ_vars.raw_file);

SBJ_vars.recon.surf_l     = [SBJ_vars.dirs.recon 'Surfaces/' SBJ_vars.SBJ '_cortex_lh.mat'];
SBJ_vars.recon.surf_r     = [SBJ_vars.dirs.recon 'Surfaces/' SBJ_vars.SBJ '_cortex_rh.mat'];
SBJ_vars.recon.elec_pat   = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_acpc_f.mat'];
SBJ_vars.recon.elec_mni_v = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_mni_frv.mat'];
SBJ_vars.recon.elec_mni_s = [];
SBJ_vars.recon.fs_T1      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_T1.mgz'];
SBJ_vars.recon.fs_DK      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_aparc+aseg.mgz'];
SBJ_vars.recon.fs_Dx      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_aparc.a2009s+aseg.mgz'];

%--------------------------------------
% Channel Selection
%--------------------------------------
SBJ_vars.ch_lab.probes     = {'LAM','LHH','LTH','LIN','LOF','LAC','LPRE',...
                              'RAM','RHH','RTH','ROF','RIN','RAC','RPRE'};
SBJ_vars.ch_lab.probe_type = {'seeg','seeg','seeg','seeg','seeg','seeg','seeg',...
                              'seeg','seeg','seeg','seeg','seeg','seeg','seeg'};
SBJ_vars.ch_lab.ref_type   = {'BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP'};
if ~all(numel(SBJ_vars.ch_lab.probes)==[numel(SBJ_vars.ch_lab.probe_type) numel(SBJ_vars.ch_lab.ref_type)]);
    error('probes ~= type+ref');end;
SBJ_vars.ch_lab.ROI        = {'all'};%'LAC*','LIN*','LOF*','LPRE*','RAC*','RIN*','ROF*','RPRE*'};
SBJ_vars.ch_lab.eeg_ROI    = {'CZ'};

%SBJ_vars.ch_lab.prefix = ''; %'POL ';    % before every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.suffix = ''; %'-Ref';    % after every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.mislabel = {};% {{'ch1bad','ch1good'}};

SBJ_vars.ch_lab.ref_exclude = {}; %exclude from the CAR
SBJ_vars.ch_lab.bad = {...
    'LTH1','LTH2','LHH1','LHH2','LHH3','LHH4','LAM1','LAM2','LAM3','LHH5',...% epileptic
    'RTH2','RTH3','RHH2','RHH3','RHH4','RHH5','RAM1','RAM2','RAM3',...% epileptic
    'LTH3','LHH6','LAM4','RAM4','RHH1','RTH1',...% spread (addeed 3/2/19 viewing)
    'RTH7','RTH8',...% spiking artifact
    'RTH10','RIN9','RIN10','RPRE10',...% out of brain (based on my eye, maybe add RAC10 and LIN10 since they're edge)
    'DC03','DC04','G','EKG','REF','E'...% not real data
    };
% bad_codes: 1 = toss (epileptic or bad); 2 = suspicious; 3 = out of brain; 0 = junk
SBJ_vars.ch_lab.bad_type = {'bad','sus','out'};
SBJ_vars.ch_lab.bad_code = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 0 0 0 0 0 0];
if numel(SBJ_vars.ch_lab.bad)~=numel(SBJ_vars.ch_lab.bad_code);error('bad ~= bad_code');end
SBJ_vars.ch_lab.eeg = {'FPZ','CZ','OZ','C3','C4'};
% SBJ_vars.ch_lab.CZ_lap_ref = {};
SBJ_vars.ch_lab.eog = {'LUC','LLC','RUC','RLC'};
SBJ_vars.ch_lab.photod = {'DC01'};

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
SBJ_vars.analysis_time = {};
SBJ_vars.ignore_trials = {[]};
if numel(SBJ_vars.analysis_time) ~= numel(SBJ_vars.raw_file) || ...
        numel(SBJ_vars.raw_file) ~= numel(SBJ_vars.block_name) || ...
        numel(SBJ_vars.ignore_trials) ~= numel(SBJ_vars.raw_file)
    error('Mismatch number of runs to concatenate!');
end

%--------------------------------------
% Artifact Rejection Parameters
%--------------------------------------
% SBJ_vars.artifact_params.std_limit_raw = 7;
% SBJ_vars.artifact_params.hard_threshold_raw = 800;

% SBJ_vars.artifact_params.std_limit_diff = 7;
% SBJ_vars.artifact_params.hard_threshold_diff = 100;

%--------------------------------------
% Trials to Reject
%--------------------------------------
% SBJ_vars.trial_reject_n = [];
