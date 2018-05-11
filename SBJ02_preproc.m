function SBJ02_preproc(SBJ,pipeline_id)
%% Preprocess data using fieldtrip
% Inputs:
%   SBJ [str]- name of the subject
%   pipeline_id [str] - name of the pipeline containing proc_vars struct

% Parameters
psd_bp       = 'yes';          % plot psds after filtering?
psd_reref    = 'yes';          % plot psds after rereferencing?
psd_line     = 'yes';          % plot psds after filtering out line noise?
psd_fig_type = 'jpg';

% Directories
addpath('/home/knight/hoycw/PRJ_Error/scripts/');
addpath('/home/knight/hoycw/PRJ_Error/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Load data
fprintf('============== Loading Data %s ==============\n',SBJ);
eval(['run /home/knight/hoycw/PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Error/scripts/proc_vars/' pipeline_id '_proc_vars.m']);

import_filename = [SBJ_vars.dirs.import SBJ '_' num2str(proc_vars.resample_freq) 'hz.mat'];
load(import_filename);
data_orig = data;

%% High and low pass data
fprintf('============== High and Low Pass Filtering %s ==============\n',SBJ);
cfg           = [];
cfg.demean    = proc_vars.demean_yn;
cfg.lpfilter  = proc_vars.lp_yn;
cfg.lpfreq    = proc_vars.lp_freq;
cfg.hpfilter  = proc_vars.hp_yn;
cfg.hpfreq    = proc_vars.hp_freq;
if isfield(proc_vars,'hp_order')
    cfg.hpfiltord = proc_vars.hp_order;
end
data = ft_preprocessing(cfg,data);
data_bp = data;

if strcmp(psd_bp,'yes')
    psd_dir = strcat(SBJ_vars.dirs.preproc,'PSDs/bp/');
    if ~exist(psd_dir,'dir')
        mkdir(psd_dir);
    end
    fn_plot_PSD_1by1_compare_save(data_orig.trial{1},data.trial{1},data_orig.label,data.label,...
        data.fsample,strcat(psd_dir,SBJ,'_PSD_bp'),'orig','bp',psd_fig_type);
end
clear data_orig

%% Rerefrence
fprintf('============== Re-Referencing %s ==============\n',SBJ);
if numel(SBJ_vars.ref_types)~=numel(SBJ_vars.ch_lab.probes)
    error('ERROR: Mismatched number of probes and reference types in SBJ_vars');
end
left_out_ch = {};
for d = 1:numel(SBJ_vars.ch_lab.probes)
    cfg = [];
    cfg.channel = ft_channelselection(strcat(SBJ_vars.ch_lab.probes{d},'*'), data.label);
    probe_data = ft_selectdata(cfg,data);   % Grab data from this probe to plot in PSD comparison
    
    % Create referencing scheme
    if strcmp(SBJ_vars.ref_types{d},'BP')
        cfg.montage.labelold = cfg.channel;
        [cfg.montage.labelnew, cfg.montage.tra, left_out_ch{d}] = fn_create_ref_scheme_bipolar(cfg.channel);
    elseif strcmp(SBJ_vars.ref_types{d},'CAR')
        cfg.reref      = 'yes';
        cfg.refchannel = setdiff(probe_data.label,SBJ_vars.ref_exclude);
        cfg.refmethod  = 'avg';
        left_out_ch{d} = {};    % CAR is applied to all channels
    else
        error(strcat('ERROR: Unrecognized reference type ',SBJ_vars.ref_types{d},...
            ' for probe ',SBJ_vars.ch_lab.probes{d}));
    end
    data_reref{d} = ft_preprocessing(cfg, data);
    
    if strcmp(psd_reref,'yes')
        psd_dir = strcat(SBJ_vars.dirs.preproc,'PSDs/bp.reref/');
        if ~exist(psd_dir,'dir')
            mkdir(psd_dir);
        end
        fn_plot_PSD_1by1_compare_save(probe_data.trial{1},data_reref{d}.trial{1},...
            probe_data.label,data_reref{d}.label,data_reref{d}.fsample,...
            strcat(psd_dir,SBJ,'_PSD_bp.reref'),'bp','bp.reref',psd_fig_type);
    end
end
clear data_bp probe_data

%Concatenate together again
cfg = [];
% cfg.appendsens = 'yes';
data = ft_appenddata(cfg,data_reref{:});
% Somehow, data.fsample is lost in certain cases here (new ft version I think):
data.fsample = data_reref{1}.fsample;
data_reref = data;

% Print left out channels, add to SBJ_vars
if ~isempty([left_out_ch{:}])
    fprintf('=============================================================\n');
    fprintf('WARNING: %i channels left out!\n',numel([left_out_ch{:}]));
    left_out_ch{:}
    fprintf('Consider adding these to SBJ_vars!\n');
    fprintf('=============================================================\n');
end
% SBJ_vars.ch_lab.left_out = [left_out_ch{:}];

%% Filter out line noise
fprintf('============== Filtering Line Noise %s via %s ==============\n',SBJ,proc_vars.notch_type);

if strcmp(proc_vars.notch_type,'dft')
    cfg           = [];
    cfg.dftfilter = 'yes'; % line noise removal using discrete fourier transform
    cfg.dftfreq   = SBJ_vars.notch_freqs;
    data = ft_preprocessing(cfg,data);
elseif strcmp(proc_vars.notch_type,'bandstop')
    % Calculate frequency limits
    bs_freq_lim = NaN([length(SBJ_vars.notch_freqs) 2]);
    for f_ix = 1:length(SBJ_vars.notch_freqs)
        bs_freq_lim(f_ix,:) = fn_freq_lim_from_CFBW(SBJ_vars.notch_freqs(f_ix),SBJ_vars.bs_width);
    end
    cfg          = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = bs_freq_lim;
    data = ft_preprocessing(cfg,data);
elseif strcmp(proc_vars.notch_type,'cleanline')
    data = fn_cleanline(data,SBJ_vars.notch_freqs);
else
    error('ERROR: proc_vars.notch_type type not in [dft, bandstop, cleanline]');
end

if strcmp(psd_line,'yes')
    psd_dir = strcat(SBJ_vars.dirs.preproc,'PSDs/',pipeline_id,'/');
    if ~exist(psd_dir,'dir')
        mkdir(psd_dir);
    end
    fn_plot_PSD_1by1_compare_save(data_reref.trial{1},data.trial{1},data_reref.label,data.label,...
        data_reref.fsample,strcat(psd_dir,SBJ,'_PSD_',pipeline_id),'bp.reref',pipeline_id,psd_fig_type);
end
clear data_reref

%% Save data
output_filename = strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',pipeline_id,'.mat');
fprintf('============== Saving %s ==============\n',output_filename);
save(output_filename, '-v7.3', 'data');

end
