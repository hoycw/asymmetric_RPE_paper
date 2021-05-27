%% IR94 Processing Variables
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
if isempty(strfind(path,'fieldtrip')); addpath(ft_dir); ft_defaults; end

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'IR94';
SBJ_vars.raw_file = {'IR94_error_raw_R1.mat','IR94_error_raw_R2.mat'};
SBJ_vars.block_name = {'R1','R2'};
SBJ_vars.low_srate  = [500,500];

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
SBJ_vars.dirs.nlx     = {[SBJ_vars.dirs.raw 'TT_run1_Hookup_2020-02-18_20-09-27/'],...
                         [SBJ_vars.dirs.raw 'TT_run2_SunPM_2020-02-23_15-07-11/']};

SBJ_vars.dirs.raw_filename = strcat(SBJ_vars.dirs.raw,SBJ_vars.raw_file);

SBJ_vars.recon.surf_l     = [SBJ_vars.dirs.recon 'Surfaces/' SBJ_vars.SBJ '_cortex_lh.mat'];
SBJ_vars.recon.surf_r     = [SBJ_vars.dirs.recon 'Surfaces/' SBJ_vars.SBJ '_cortex_rh.mat'];
SBJ_vars.recon.wm_l       = [SBJ_vars.dirs.recon 'Surfaces/' SBJ_vars.SBJ '_wm_lh.mat'];
SBJ_vars.recon.wm_r       = [SBJ_vars.dirs.recon 'Surfaces/' SBJ_vars.SBJ '_wm_rh.mat'];
SBJ_vars.recon.infl_l     = [SBJ_vars.dirs.recon 'Surfaces/' SBJ_vars.SBJ '_infl_lh.mat'];
SBJ_vars.recon.infl_r     = [SBJ_vars.dirs.recon 'Surfaces/' SBJ_vars.SBJ '_infl_rh.mat'];
SBJ_vars.recon.elec_pat   = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_acpc_....mat'];
SBJ_vars.recon.elec_mni_v = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_mni_frv.mat'];
SBJ_vars.recon.elec_mni_s = [];%[SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_mni_s.mat'];
SBJ_vars.recon.fs_T1      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_T1.mgz'];
SBJ_vars.recon.fs_DK      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_aparc+aseg.mgz'];
SBJ_vars.recon.fs_Dx      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_aparc.a2009s+aseg.mgz'];

%--------------------------------------
% Channel Selection
%--------------------------------------
SBJ_vars.ch_lab.probes     = {'RAM','RHH','RTH','RPT','RPC'...
                              'PCU','RAI','RPI','RBT','CUN'};
SBJ_vars.ch_lab.probe_type = {'seeg','seeg','seeg','seeg','seeg',...
                              'seeg','seeg','seeg','seeg','seeg'}; % 'ecog' or 'seeg'
SBJ_vars.ch_lab.ref_type   = {'BP','BP','BP','BP','BP','BP','BP','BP','BP','BP'}; % 'BP', 'CAR', 'CARall', 'CMR'
if ~all(numel(SBJ_vars.ch_lab.probes)==[numel(SBJ_vars.ch_lab.probe_type) numel(SBJ_vars.ch_lab.ref_type)]); error('probes ~= type+ref');end;
SBJ_vars.ch_lab.nlx        = [1,1,1,1,1,1,1,0,0,0]; % 0/1 if these macros are also in NLX system
SBJ_vars.ch_lab.ROI        = {'all'}; % 'all' or subset of probes
SBJ_vars.ch_lab.eeg_ROI    = {}; % 'all' or subset of channels
SBJ_vars.ch_lab.wires      = {'mram','mrhh','mrth','mpcu','mrpc','mrpt'}; % list of microwire bundles {'mram','mlac',...}
SBJ_vars.ch_lab.wire_type  = {'su','su','su','su','su','su'}; % 'su' (anything else?)
SBJ_vars.ch_lab.wire_ref   = {'','','','','',''}; % theoretically which one is the reference, but unused?
SBJ_vars.ch_lab.wire_ROI   = {'all'}; % which microwires to use

%SBJ_vars.ch_lab.prefix = 'POL ';    % before every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.suffix = '-Ref';    % after every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.mislabel = {{'RLT12','FPG12'},{'IH;L8','IHL8'}};

SBJ_vars.ch_lab.nlx_suffix   = {'_0006',''}; % session suffix for files, e.g., '_0007'
SBJ_vars.ch_lab.nlx_nk_align = {'RAM4','RAM5'}; % shared macros fro clinical-NLX alignment, 1 for unipolar or 2 for bipolar
SBJ_vars.nlx_macro_inverted  = [1,1]; % 0/1 is NLX recorded inverted? usually 1
SBJ_vars.nlx_analysis_time   = {{[0 1640]},{[0 1371.5]}}; % cut to time in NLX photodiode (e.g., discontinuities)

SBJ_vars.ch_lab.ref_exclude = {}; %exclude from the CAR
SBJ_vars.ch_lab.bad = {...
    'RBT1','RHH6',... % loose channel
    'RHH7','RHH8','RTH8','RAI7','RAI8','RAI9','RPI9','RAI10','RPI10','RBT10',... % out of brain maybe RTH7
    'EKG',... % EKG channel
    'DC01','DC02','DC03','DC04',... % empty analog channels
    'Events','GND','XREF','Mark1','Mark2','E'... % junk channels
    };
% I suspect RAM1 is a source that spreads to RAM2, RHH1, and RTH1
% bad_codes: 1 = toss (epileptic or bad); 2 = suspicious; 3 = out of brain; 0 = junk
SBJ_vars.ch_lab.bad_type = {'bad','sus','out'};
SBJ_vars.ch_lab.bad_code = [1,1,3,3,3,3,3,3,3,3,3,3,0,0,0,0,0,0,0,0,0,0,0];
if numel(SBJ_vars.ch_lab.bad)~=numel(SBJ_vars.ch_lab.bad_code);error('bad ~= bad_code');end
SBJ_vars.ch_lab.eeg = {'FZ','CZ','OZ','C3','C4'}; % scalp channel labels
% SBJ_vars.ch_lab.CZ_lap_ref = {}; % reference channels for scalp laplacian
SBJ_vars.ch_lab.eog = {'LUC','LLC','RLC','RUC'}; % EOG channel labels
SBJ_vars.ch_lab.photod = {'photo1'}; % photodiode label
SBJ_vars.photo_inverted = [1,1]; % 0/1
SBJ_vars.ch_lab.mic    = {'mic1'}; % microphone label

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
SBJ_vars.analysis_time = {{[0 1700]},{[0 1440]}};
SBJ_vars.ignore_trials = {[],[]};
if numel(SBJ_vars.analysis_time) ~= numel(SBJ_vars.raw_file) || ...
        numel(SBJ_vars.raw_file) ~= numel(SBJ_vars.block_name) || ...
        numel(SBJ_vars.ignore_trials) ~= numel(SBJ_vars.raw_file)
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
