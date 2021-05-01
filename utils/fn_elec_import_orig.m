function fn_elec_import_orig(SBJ,proc_id,view_space,reg_type,reref)
%% Convert recon pipeline elec struct to personal format
%   Fix labels, select imported channels (ROI, cleaning), sort alphanumerically
%   Rereference: if 1, does bipolar (no difference for ECoG)
%       Goes from smallest to largest (ELEC1-ELEC2, ELEC2-ELEC3, etc.)
%   adds channel type, hemisphere
% INPUTS:
%   SBJ [str] - name of subject
%   proc_id [str] - name of analysis pipeline
%   view_space [str] - {'pat','mni'} select patient native or mni group space
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   reref [0/1] - apply preprocessing re-referencing?

% Set up paths
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load variables
eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/proc_vars/' proc_id '_vars.m']);

%% Load Elec struct
[elec] = fn_load_elec_orig(SBJ,view_space,reg_type);

% Create correct extension for saving
elec_ext = view_space;
if ~isempty(reg_type)
    if strcmp(view_space,'pat')
        error('pat view_space should have empty reg_type');
    end
    elec_ext = [elec_ext '_' reg_type];
end
if ~reref
    elec_ext = [elec_ext '_orig'];
end

%% Channel Label Corrections
% Strip Pre/Suffix if Necessary
for ch_ix = 1:numel(elec.label)
    if isfield(SBJ_vars.ch_lab,'prefix')
        elec.label{ch_ix} = strrep(elec.label{ch_ix},SBJ_vars.ch_lab.prefix,'');
    end
    if isfield(SBJ_vars.ch_lab,'suffix')
        elec.label{ch_ix} = strrep(elec.label{ch_ix},SBJ_vars.ch_lab.suffix,'');
    end
end
% Fix any mislabeled channels
if isfield(SBJ_vars.ch_lab,'mislabel')
    for ch_ix = 1:numel(SBJ_vars.ch_lab.mislabel)
        % Future edit: search for the bad label across data, eeg, evnt
        elec.label(strcmp(elec.label,SBJ_vars.ch_lab.mislabel{ch_ix}(1))) = SBJ_vars.ch_lab.mislabel{ch_ix}(2);
    end
end

%% Load data
if numel(SBJ_vars.raw_file)>1
    block_suffix = strcat('_',SBJ_vars.block_name{1});
else
    block_suffix = SBJ_vars.block_name{1};   % should just be ''
end
if any(SBJ_vars.low_srate)
    import_fname = [SBJ_vars.dirs.import SBJ '_',num2str(SBJ_vars.low_srate(1)),'hz',block_suffix,'.mat'];
else
    import_fname = [SBJ_vars.dirs.import SBJ '_',num2str(proc.resample_freq),'hz',block_suffix,'.mat'];
end
load(import_fname);

% % Original (single electrode) labels
% import  = load([SBJ_vars.dirs.import SBJ '_1000hz.mat']);
% raw_lab = import.data.label;

%% Select imported channels
cfg = []; cfg.channel = data.label;
elec = fn_select_elec(cfg,elec);

% Order them to match data.label
elec = fn_reorder_elec(elec, data.label);
SBJ_vars.ch_lab.probes = sort(SBJ_vars.ch_lab.probes);  %alphabetical, like preproc

%% Apply montage per probe
left_out_ch = {};
elec_labels = {};
elec_types  = {};
danger_name = false([1 numel(SBJ_vars.ch_lab.probes)]);
name_holder = cell([2 numel(SBJ_vars.ch_lab.probes)]);
elec_reref  = cell([1 numel(SBJ_vars.ch_lab.probes)]);
for d = 1:numel(SBJ_vars.ch_lab.probes)
    cfg = [];
    cfg.channel = ft_channelselection(strcat(SBJ_vars.ch_lab.probes{d},'*'), elec.label);
    probe_elec  = fn_select_elec(cfg,elec);
    %     probe_data = ft_selectdata(cfg,data);   % Grab data from this probe to plot in PSD comparison
    %     probe_data.elec = fn_elec_ch_select(elec,cfg.channel);
    
    % Check if the names of these elecs will cause problems
    eeg1010_match = strfind(probe_elec.label,'AF');
    if ~isempty([eeg1010_match{:}])
        danger_name(d)   = true;
        name_holder{1,d} = probe_elec.label;
        name_holder{2,d} = fn_generate_random_strings(numel(probe_elec.label),'',10);
        probe_elec.label = name_holder{2,d};
    end
    
    % Create referencing scheme
    if reref && strcmp(SBJ_vars.ch_lab.ref_type{d},'BP')
        cfg.montage.labelold = cfg.channel;
        [cfg.montage.labelnew, cfg.montage.tra, left_out_ch{d}] = fn_create_ref_scheme_bipolar(cfg.channel);
        cfg.updatesens = 'yes';
        elec_reref{d} = ft_apply_montage(probe_elec, cfg.montage);%, 'feedback', 'none', 'keepunused', 'no', 'balancename', bname);
        %     data_reref{d} = ft_preprocessing(cfg, probe_data);
    else
        elec_reref{d} = probe_elec;
    end
    if d==1
        elec_labels = elec_reref{d}.label;
        elec_types  = repmat(SBJ_vars.ch_lab.probe_type(d),size(elec_reref{d}.label));
    else
        elec_labels = cat(find(size(elec_labels)>1), elec_labels, elec_reref{d}.label);
        elec_types  = cat(find(size(elec_types)>1), elec_types, repmat(SBJ_vars.ch_lab.probe_type(d),size(elec_reref{d}.label)));
    end
end

% Recombine
cfg = [];
elec = ft_appendsens(cfg,elec_reref{:});

% Re-label any problematic channel labels
if any(danger_name)
    for d_ix = find(danger_name)
        for s_ix = 1:numel(name_holder{2,d_ix})
            elec.label{strcmp(elec.label,name_holder{2,2}{s_ix})} = name_holder{1,d_ix}{s_ix};
        end
    end
end

% Remove non-ROI (bad) electrodes (shouldn't do anything if non-reref)
cfgs = [];
cfgs.channel = SBJ_vars.ch_lab.ROI;
elec = fn_select_elec(cfgs,elec);

% Sort Alphanumerically
elec = fn_reorder_elec(elec,'');

%% Add Channel Types
elec.type = 'ieeg';
for e = 1:numel(elec.chantype)
    elec.chantype{e} = elec_types{strcmp(elec.label{e},elec_labels)};
end

%% Add in L/R Hemisphere
elec.hemi = repmat({'r'},size(elec.label));
for e = 1:numel(elec.label)
    if strfind(elec.label{e},'L')
        elec.hemi{e} = 'l';
    end
end

% Fix Exceptions
if strcmp(SBJ,'IR68')                                       % IR68- All probes in L
    elec.hemi = repmat({'l'},size(elec.label));
elseif strcmp(SBJ,'CP24')                                   % CP24- RLF has L in it
    rlf_ix = ~cellfun(@isempty,strfind(elec.label,'RLF'));
    elec.hemi(rlf_ix) = repmat({'r'},size(elec.label(rlf_ix)));
end

%% Save data
% Check if elec.cfg.previosu got ridiculously large, and keep only first
var_stats = whos('elec');
if var_stats.bytes>1000000
    elec.cfg = rmfield(elec.cfg,'previous');
end
output_fname = strcat(SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_',elec_ext,'.mat');
fprintf('============== Saving %s ==============\n',output_fname);
save(output_fname, '-v7.3', 'elec');

end
