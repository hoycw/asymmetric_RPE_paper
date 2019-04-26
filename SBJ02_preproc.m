function SBJ02_preproc(SBJ,proc_id)
%% Preprocess data using fieldtrip
% Inputs:
%   SBJ [str]- name of the subject
%   proc_id [str] - name of the pipeline containing proc_vars struct

% Parameters
psd_bp       = 'yes';          % plot psds after filtering?
psd_reref    = 'yes';          % plot psds after rereferencing?
psd_line     = 'yes';          % plot psds after filtering out line noise?
psd_fig_type = 'jpg';

% Directories
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath(genpath([root_dir 'PRJ_Error/scripts/utils/']));
addpath(ft_dir);
ft_defaults

%% Load data
fprintf('============== Loading Data %s ==============\n',SBJ);
eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/proc_vars/' proc_id '_proc_vars.m']);

data_all = {};
for b_ix = 1:numel(SBJ_vars.block_name)
    if numel(SBJ_vars.raw_file)>1
        block_suffix = strcat('_',SBJ_vars.block_name{b_ix});
    else
        block_suffix = SBJ_vars.block_name{b_ix};   % should just be ''
    end
    if any(SBJ_vars.low_srate)
        import_fname = [SBJ_vars.dirs.import SBJ '_',num2str(SBJ_vars.low_srate(b_ix)),'hz',block_suffix,'.mat'];
    else
        import_fname = [SBJ_vars.dirs.import SBJ '_',num2str(proc_vars.resample_freq),'hz',block_suffix,'.mat'];
    end
    load(import_fname);
    data_orig = data;
    
    %% High and low pass data
    fprintf('============== High and Low Pass Filtering %s ==============\n',SBJ);
    cfg           = [];
    cfg.demean    = proc_vars.demean_yn;
    if numel(SBJ_vars.block_name)>1
        cfg.demean = 'yes';
    end
    if proc_vars.lp_freq > data.fsample/2
        cfg.lpfilter = 'no';
    else
        cfg.lpfilter = proc_vars.lp_yn;
        cfg.lpfreq   = proc_vars.lp_freq;
    end
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
            data.fsample,strcat(psd_dir,SBJ,'_PSD_bp',block_suffix),'orig','bp',psd_fig_type);
    end
    clear data_orig
    
    %% Rereference
    fprintf('============== Re-Referencing %s ==============\n',SBJ);
    left_out_ch = {};
    if all(strcmp(SBJ_vars.ch_lab.ref_type,'CARall'))
        % CAR across all channels! (likely for CPMC data)
        cfg = [];
        cfg.reref      = 'yes';
        cfg.refchannel = setdiff(data.label,SBJ_vars.ch_lab.ref_exclude);
        cfg.refmethod  = 'avg';
        cfg.updatesens = 'yes';
        data_reref = ft_preprocessing(cfg, data);
        data = data_reref;
        % left_out_ch stays empty, CAR is applied to all channels
        
        if strcmp(psd_reref,'yes')
            psd_dir = strcat(SBJ_vars.dirs.preproc,'PSDs/bp.reref/');
            if ~exist(psd_dir,'dir')
                mkdir(psd_dir);
            end
            fn_plot_PSD_1by1_compare_save(data.trial{1},data_reref.trial{1},...
                data.label,data_reref.label,data_reref.fsample,...
                strcat(psd_dir,SBJ,'_PSD_bp.reref',block_suffix),'bp','bp.reref',psd_fig_type);
        end
    elseif numel(SBJ_vars.ch_lab.ref_type)==numel(SBJ_vars.ch_lab.probes)
        SBJ_vars.ch_lab.probes = sort(SBJ_vars.ch_lab.probes);  % keep alphabetical order
        for d = 1:numel(SBJ_vars.ch_lab.probes)
            cfg = [];
            cfg.channel = ft_channelselection(strcat(SBJ_vars.ch_lab.probes{d},'*'), data.label);
            probe_data = ft_selectdata(cfg,data);   % Grab data from this probe to plot in PSD comparison
            
            % Create referencing scheme
            if strcmp(SBJ_vars.ch_lab.ref_type{d},'BP')
                cfg.montage.labelold = cfg.channel;
                [cfg.montage.labelnew, cfg.montage.tra, left_out_ch{d}] = fn_create_ref_scheme_bipolar(cfg.channel);
            elseif any(strcmp(SBJ_vars.ch_lab.ref_type{d},{'CAR','CMR'}))
                left_out_ch{d} = {};    % CAR/CMR is applied to all channels
                cfg.reref      = 'yes';
                cfg.refchannel = setdiff(probe_data.label,SBJ_vars.ch_lab.ref_exclude);
                if strcmp(SBJ_vars.ch_lab.ref_type{d},'CMR')
                    cfg.refmethod  = 'median';
                else
                    cfg.refmethod  = 'avg';
                end
            else
                error(strcat('ERROR: Unrecognized reference type ',SBJ_vars.ch_lab.ref_type{d},...
                    ' for probe ',SBJ_vars.ch_lab.probes{d}));
            end
            cfg.updatesens = 'yes';
            data_reref{d} = ft_preprocessing(cfg, data);
            
            if strcmp(psd_reref,'yes')
                psd_dir = strcat(SBJ_vars.dirs.preproc,'PSDs/bp.reref/');
                if ~exist(psd_dir,'dir')
                    mkdir(psd_dir);
                end
                fn_plot_PSD_1by1_compare_save(probe_data.trial{1},data_reref{d}.trial{1},...
                    probe_data.label,data_reref{d}.label,data_reref{d}.fsample,...
                    strcat(psd_dir,SBJ,'_PSD_bp.reref',block_suffix),'bp','bp.reref',psd_fig_type);
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
    else
        error('ERROR: Mismatched number of probes and reference types in SBJ_vars');
    end
    
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
        cfg.dftfreq(cfg.dftfreq > data.fsample/2) = [];
        data = ft_preprocessing(cfg,data);
    elseif strcmp(proc_vars.notch_type,'bandstop')
        % Calculate frequency limits
        bs_freq_lim = NaN([length(SBJ_vars.notch_freqs) 2]);
        for f_ix = 1:length(SBJ_vars.notch_freqs)
            bs_freq_lim(f_ix,:) = fn_freq_lim_from_CFBW(SBJ_vars.notch_freqs(f_ix),SBJ_vars.bs_width);
        end
        bs_freq_lim(bs_freq_lim(:,2) > data.fsample/2, :) = [];
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
        psd_dir = strcat(SBJ_vars.dirs.preproc,'PSDs/',proc_id,'/');
        if ~exist(psd_dir,'dir')
            mkdir(psd_dir);
        end
        fn_plot_PSD_1by1_compare_save(data_reref.trial{1},data.trial{1},data_reref.label,data.label,...
            data_reref.fsample,strcat(psd_dir,SBJ,'_PSD_',proc_id,block_suffix),...
            'bp.reref',proc_id,psd_fig_type);
    end
    
    data_all{b_ix} = data;
    clear data_reref
end

%% Concatenate blocks
data = fn_concat_blocks(data_all);

%% Save data
output_filename = strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat');
fprintf('============== Saving %s ==============\n',output_filename);
save(output_filename, '-v7.3', 'data');

end
