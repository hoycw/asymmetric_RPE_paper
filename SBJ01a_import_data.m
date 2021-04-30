function SBJ01a_import_data(SBJ,proc_id)
%% SBJ01a_import_data.m
% Extract data with fieldtrip and save out by data type

if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load and preprocess the data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/proc_vars/' proc_id '_proc_vars.m'];
eval(proc_vars_cmd);

%% Process channel labels
if isfield(SBJ_vars.ch_lab,'prefix')
    for bad_ix = 1:numel(SBJ_vars.ch_lab.bad)
        SBJ_vars.ch_lab.bad{bad_ix} = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.bad{bad_ix}];
    end
    for eeg_ix = 1:numel(SBJ_vars.ch_lab.eeg)
        SBJ_vars.ch_lab.eeg{eeg_ix} = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.eeg{eeg_ix}];
    end
    for eog_ix = 1:numel(SBJ_vars.ch_lab.eog)
        SBJ_vars.ch_lab.eog{eog_ix} = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.eog{eog_ix}];
    end
    SBJ_vars.ch_lab.photod = {[SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.photod{1}]};
end
if isfield(SBJ_vars.ch_lab,'suffix')
    for bad_ix = 1:numel(SBJ_vars.ch_lab.bad)
        SBJ_vars.ch_lab.bad{bad_ix} = [SBJ_vars.ch_lab.bad{bad_ix} SBJ_vars.ch_lab.suffix];
    end
    for eeg_ix = 1:numel(SBJ_vars.ch_lab.eeg)
        SBJ_vars.ch_lab.eeg{eeg_ix} = [SBJ_vars.ch_lab.eeg{eeg_ix} SBJ_vars.ch_lab.suffix];
    end
    for eog_ix = 1:numel(SBJ_vars.ch_lab.eog)
        SBJ_vars.ch_lab.eog{eog_ix} = [SBJ_vars.ch_lab.eog{eog_ix} SBJ_vars.ch_lab.suffix];
    end
    SBJ_vars.ch_lab.photod = {[SBJ_vars.ch_lab.photod{1} SBJ_vars.ch_lab.suffix]};
end
bad_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.bad);
eeg_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.eeg);
eog_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.eog);
photod_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.photod);

for b_ix = 1:numel(SBJ_vars.block_name)
    %% Load Data
    if numel(SBJ_vars.raw_file)>1
        block_suffix = strcat('_',SBJ_vars.block_name{b_ix});
    else
        block_suffix = SBJ_vars.block_name{b_ix};   % should just be ''
    end
    if strcmp(SBJ_vars.raw_file{b_ix}(end-2:end),'mat')
        orig = load(SBJ_vars.dirs.raw_filename{b_ix});
        % Load Neural Data
        cfg = [];
        cfg.channel = [{'all'},bad_ch_neg,eeg_ch_neg,eog_ch_neg,photod_ch_neg];
        data = ft_selectdata(cfg,orig.data);
        % EEG data
        if ~isempty(SBJ_vars.ch_lab.eeg)
            cfg.channel = SBJ_vars.ch_lab.eeg;
            eeg = ft_selectdata(cfg,orig.data);
        end
        % EOG data
        if ~isempty(SBJ_vars.ch_lab.eog)
            cfg.channel = SBJ_vars.ch_lab.eog;
            eog = ft_selectdata(cfg,orig.data);
        end
        % Load event data
        if ~isfield(SBJ_vars.dirs,'nlx')
            cfg.channel = SBJ_vars.ch_lab.photod;
            evnt = ft_selectdata(cfg,orig.data);
        end
    else
        % Load Neural Data
        cfg            = [];
        cfg.dataset    = SBJ_vars.dirs.raw_filename{b_ix};
        cfg.continuous = 'yes';
        cfg.channel    = [{'all'},bad_ch_neg,eeg_ch_neg,eog_ch_neg,photod_ch_neg];
        data = ft_preprocessing(cfg);
        
        % Load EEG data
        if ~isempty(SBJ_vars.ch_lab.eeg)
            cfg.channel = SBJ_vars.ch_lab.eeg;
            eeg = ft_preprocessing(cfg);
        end
        
        % Load EOG data
        if ~isempty(SBJ_vars.ch_lab.eog)
            cfg.channel = SBJ_vars.ch_lab.eog;
            eog = ft_preprocessing(cfg);
        end
        
        % Load event data
        if ~isfield(SBJ_vars.dirs,'nlx')
            cfg.channel = SBJ_vars.ch_lab.photod;
            evnt = ft_preprocessing(cfg);
        end
    end
    % Toss 'EDF Annotations' Channel (if it exists)
    if any(strcmp(data.label,'EDF Annotations'))
        cfg_s = [];
        cfg_s.channel = {'all','-EDF Annotations'};
        data = ft_selectdata(cfg_s,data);
    end
    
    %% Cut to analysis_time
    if ~isempty(SBJ_vars.analysis_time{b_ix})
        for epoch_ix = 1:length(SBJ_vars.analysis_time{b_ix})
            epoch_len{epoch_ix} = diff(SBJ_vars.analysis_time{b_ix}{epoch_ix});
            cfg_trim = [];
            cfg_trim.latency = SBJ_vars.analysis_time{b_ix}{epoch_ix};
            
            data_pieces{epoch_ix} = ft_selectdata(cfg_trim, data);
            if ~isfield(SBJ_vars.dirs,'nlx')
                evnt_pieces{epoch_ix} = ft_selectdata(cfg_trim, evnt);
            end
            if ~isempty(SBJ_vars.ch_lab.eeg)
                eeg_pieces{epoch_ix}  = ft_selectdata(cfg_trim, eeg);
            end
            if ~isempty(SBJ_vars.ch_lab.eog)
                eog_pieces{epoch_ix}  = ft_selectdata(cfg_trim, eog);
            end
        end
        % Stitch them back together
        data = data_pieces{1};
        data.time{1} = data.time{1}-SBJ_vars.analysis_time{b_ix}{1}(1);
        if ~isfield(SBJ_vars.dirs,'nlx')
            evnt = evnt_pieces{1};
            evnt.time{1} = evnt.time{1}-SBJ_vars.analysis_time{b_ix}{1}(1);
        end
        if ~isempty(SBJ_vars.ch_lab.eeg)
            eeg = eeg_pieces{1};
            eeg.time{1} = eeg.time{1}-SBJ_vars.analysis_time{b_ix}{1}(1);
        end
        if ~isempty(SBJ_vars.ch_lab.eog)
            eog = eog_pieces{1};
            eog.time{1} = eog.time{1}-SBJ_vars.analysis_time{b_ix}{1}(1);
        end
        if length(SBJ_vars.analysis_time{b_ix})>1
            for epoch_ix = 2:length(SBJ_vars.analysis_time{b_ix})
                data.trial{1} = horzcat(data.trial{1},data_pieces{epoch_ix}.trial{1});
                data.time{1} = horzcat(data.time{1},data_pieces{epoch_ix}.time{1}-...
                    SBJ_vars.analysis_time{b_ix}{epoch_ix}(1)+data.time{1}(end)+data.time{1}(2));
                
                if ~isfield(SBJ_vars.dirs,'nlx')
                    evnt.trial{1} = horzcat(evnt.trial{1},evnt_pieces{epoch_ix}.trial{1});
                    evnt.time{1} = horzcat(evnt.time{1},evnt_pieces{epoch_ix}.time{1}-...
                        SBJ_vars.analysis_time{b_ix}{epoch_ix}(1)+evnt.time{1}(end)+evnt.time{1}(2));
                end
                
                if ~isempty(SBJ_vars.ch_lab.eeg)
                    eeg.trial{1} = horzcat(eeg.trial{1},eeg_pieces{epoch_ix}.trial{1});
                    eeg.time{1} = horzcat(eeg.time{1},eeg_pieces{epoch_ix}.time{1}-...
                        SBJ_vars.analysis_time{b_ix}{epoch_ix}(1)+eeg.time{1}(end)+eeg.time{1}(2));
                end
                
                if ~isempty(SBJ_vars.ch_lab.eog)
                    eog.trial{1} = horzcat(eog.trial{1},eog_pieces{epoch_ix}.trial{1});
                    eog.time{1} = horzcat(eog.time{1},eog_pieces{epoch_ix}.time{1}-...
                        SBJ_vars.analysis_time{b_ix}{epoch_ix}(1)+eog.time{1}(end)+eog.time{1}(2));
                end
            end
        end
    end
    
    %% Resample data
    if strcmp(proc_vars.resample_yn,'yes') && (data.fsample > proc_vars.resample_freq)
        cfg = [];
        cfg.resamplefs = proc_vars.resample_freq;
        cfg.detrend = 'no';
        data = ft_resampledata(cfg, data);
        if ~isfield(SBJ_vars.dirs,'nlx')
            evnt = ft_resampledata(cfg, evnt);
        end
        if ~isempty(SBJ_vars.ch_lab.eeg)
            eeg = ft_resampledata(cfg, eeg);
        end
        if ~isempty(SBJ_vars.ch_lab.eog)
            eog = ft_resampledata(cfg, eog);
        end
    end
    
    %% Final Channel Label Corrections
    % Strip Pre/Suffix if Necessary
    for ch_ix = 1:numel(data.label)
        if isfield(SBJ_vars.ch_lab,'prefix')
            data.label{ch_ix} = strrep(data.label{ch_ix},SBJ_vars.ch_lab.prefix,'');
        end
        if isfield(SBJ_vars.ch_lab,'suffix')
            data.label{ch_ix} = strrep(data.label{ch_ix},SBJ_vars.ch_lab.suffix,'');
        end
    end
    if ~isempty(SBJ_vars.ch_lab.eeg)
        for eeg_ix = 1:numel(eeg.label)
            if isfield(SBJ_vars.ch_lab,'prefix')
                eeg.label{eeg_ix} = strrep(eeg.label{eeg_ix},SBJ_vars.ch_lab.prefix,'');
            end
            if isfield(SBJ_vars.ch_lab,'suffix')
                eeg.label{eeg_ix} = strrep(eeg.label{eeg_ix},SBJ_vars.ch_lab.suffix,'');
            end
        end
    end
    if ~isempty(SBJ_vars.ch_lab.eog)
        for eog_ix = 1:numel(eog.label)
            if isfield(SBJ_vars.ch_lab,'prefix')
                eog.label{eog_ix} = strrep(eog.label{eog_ix},SBJ_vars.ch_lab.prefix,'');
            end
            if isfield(SBJ_vars.ch_lab,'suffix')
                eog.label{eog_ix} = strrep(eog.label{eog_ix},SBJ_vars.ch_lab.suffix,'');
            end
        end
    end
    if ~isfield(SBJ_vars.dirs,'nlx')
        if isfield(SBJ_vars.ch_lab,'prefix')
            evnt.label{1} = strrep(evnt.label{1},SBJ_vars.ch_lab.prefix,'');
        end
        if isfield(SBJ_vars.ch_lab,'suffix')
            evnt.label{1} = strrep(evnt.label{1},SBJ_vars.ch_lab.suffix,'');
        end
    end
    
    % Fix any mislabeled channels
    if isfield(SBJ_vars.ch_lab,'mislabel')
        for ch_ix = 1:numel(SBJ_vars.ch_lab.mislabel)
            if any(strcmp(data.label,SBJ_vars.ch_lab.mislabel{ch_ix}(1)))
                data.label(strcmp(data.label,SBJ_vars.ch_lab.mislabel{ch_ix}(1))) = SBJ_vars.ch_lab.mislabel{ch_ix}(2);
            elseif any(strcmp(eeg.label,SBJ_vars.ch_lab.mislabel{ch_ix}(1)))
                eeg.label(strcmp(eeg.label,SBJ_vars.ch_lab.mislabel{ch_ix}(1))) = SBJ_vars.ch_lab.mislabel{ch_ix}(2);
            elseif any(strcmp(eog.label,SBJ_vars.ch_lab.mislabel{ch_ix}(1)))
                eog.label(strcmp(eog.label,SBJ_vars.ch_lab.mislabel{ch_ix}(1))) = SBJ_vars.ch_lab.mislabel{ch_ix}(2);
            else
                error(['Could not find mislabeled channel: ' SBJ_vars.ch_lab.mislabel{ch_ix}(1)]);
            end
        end
    end
    
    % Check that no bad channels snuck through (e.g., mislabeled)
    for ch_ix = 1:numel(SBJ_vars.ch_lab.bad)
        if any(strcmp(SBJ_vars.ch_lab.bad{ch_ix},data.label))
            error(['ERROR: bad channel still in imported data after re-labeling: ' SBJ_vars.ch_lab.bad{ch_ix}]);
        end
    end
    
    % Sort channels alphabetically
    data = fn_reorder_data(data, {});
    if ~isempty(SBJ_vars.ch_lab.eeg)
        for ch_ix = 1:numel(SBJ_vars.ch_lab.bad)
            if any(strcmp(SBJ_vars.ch_lab.bad{ch_ix},eeg.label))
                error(['ERROR: bad channel still in imported data after re-labeling: ' SBJ_vars.ch_lab.bad{ch_ix}]);
            end
        end
        eeg = fn_reorder_data(eeg, sort(eeg.label));
    end
    if ~isempty(SBJ_vars.ch_lab.eog)
        for ch_ix = 1:numel(SBJ_vars.ch_lab.bad)
            if any(strcmp(SBJ_vars.ch_lab.bad{ch_ix},eog.label))
                error(['ERROR: bad channel still in imported data after re-labeling: ' SBJ_vars.ch_lab.bad{ch_ix}]);
            end
        end
        eog = fn_reorder_data(eog, sort(eog.label));
    end
    
    %% Save data
    nrl_out_fname = strcat(SBJ_vars.dirs.import,SBJ,'_',num2str(data.fsample),'hz',block_suffix,'.mat');
    save(nrl_out_fname, '-v7.3', 'data');
    
    if ~isempty(SBJ_vars.ch_lab.eeg)
        eeg_out_fname = strcat(SBJ_vars.dirs.import,SBJ,'_eeg_',num2str(eeg.fsample),'hz',block_suffix,'.mat');
        save(eeg_out_fname, '-v7.3', 'eeg');
    end
    
    if ~isempty(SBJ_vars.ch_lab.eog)
        eog_out_fname = strcat(SBJ_vars.dirs.import,SBJ,'_eog_',num2str(eog.fsample),'hz',block_suffix,'.mat');
        save(eog_out_fname, '-v7.3', 'eog');
    end
    
    if ~isfield(SBJ_vars.dirs,'nlx')
        evnt_out_fname = strcat(SBJ_vars.dirs.import,SBJ,'_evnt_',num2str(evnt.fsample),'hz',block_suffix,'.mat');
        save(evnt_out_fname, '-v7.3', 'evnt');
    end
    
    clear data evnt eeg eog
end

end
