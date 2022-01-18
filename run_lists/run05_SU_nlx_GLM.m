%% Preprocessing Pipeline
% This script should be run in sections. Functions/scripts with the SBJ##
% prefix can be run automatically, and all other sections should be
% manually editted for each dataset.
clear all; close all;

%% Check which root directory
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'PRJ_Error/scripts/']);
addpath([root_dir 'PRJ_Error/scripts/utils/']);
addpath(genpath('/Users/colinhoy/Code/Apps/wave_clus/'));
addpath(ft_dir);
ft_defaults

%% Step 0 - Processing Variables
% SBJ = 'IR';

proc_id = 'SU_nlx';
eval(['run ' root_dir 'PRJ_Error/scripts/proc_vars/' proc_id '_vars.m']);

%% ========================================================================
%   Step 1- Load SBJ and Processing Variable Structures
%  ========================================================================
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

%% ========================================================================
%   Step 2- Quick Import and Processing for Data Cleaning/Inspection
%  ========================================================================
SU00_extract_evnt(SBJ)
% SU00_extract_data(SBJ,proc.resample_freq??)
% block_prefix = '';
% !!! outdated command SBJ00_cleaning_prep(SBJ,SBJ_vars.raw_file,proc.plot_psd,block_prefix);

%% ========================================================================
% If you want to spike sort the .nse files:
% SU_convert_nse_wave_clus

% Spike Sorting
for b_ix = 1:numel(SBJ_vars.block_name)
    SU_wave_clus_auto(SBJ,[SBJ_vars.dirs.nlx{b_ix} 'micro/'])
end

%% ========================================================================
%   Step 3- Import Data, Resample, and Save Individual Data Types
%  ========================================================================
% FILE TOO BIG, RUNNING THIS VIA SGE
% SBJ01_import_data(SBJ,proc.resample_freq);

%% ========================================================================
%   Step 4a- Manually Clean Photodiode Trace: Load & Plot
%  ========================================================================
% Load data
for b_ix = 1:numel(SBJ_vars.block_name)
    if numel(SBJ_vars.raw_file)>1
        block_suffix = strcat('_',SBJ_vars.block_name{b_ix});
    else
        block_suffix = SBJ_vars.block_name{b_ix};   % should just be ''
    end
    evnt_filename = strcat(SBJ_vars.dirs.import,SBJ,'_evnt',block_suffix,'.mat');
    load(evnt_filename);
    [evnt, hdr] = fn_format_data_ft2KLA(evnt);
    
    % Plot event channels
    plot(linspace(0,hdr.length_in_seconds,hdr.n_samples), evnt(2,:)); hold on;
    plot(linspace(0,hdr.length_in_seconds,hdr.n_samples), evnt(1,:));
    
    % Get baseline value:
    % mean(evnt(1,start*hdr.sample_rate:stop*hdr.sample_rate))
    
    %% ========================================================================
    %   Step 4b- Manually Clean Photodiode Trace: Mark Sections to Correct
    %  ========================================================================
    % Create correction times and values in a separate file in ~/PRJ_Error/scripts/SBJ_evnt_clean/
    SBJ_evnt_clean_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_evnt_clean/' SBJ '_evnt_clean_params',block_suffix,'.m'];
    eval(SBJ_evnt_clean_cmd);
    
    %% ========================================================================
    %   Step 4c- Manually Clean Photodiode Trace: Apply Corrections
    %  ========================================================================
    photod_ix = strmatch('photo',hdr.channel_labels);
    mic_ix = strmatch('mic',hdr.channel_labels);
    
    % If contionuous drift, use cftool to find slop on empty period and subtract
    % drift_model = slope*[1:size(evnt,2)];
    % evnt(photod_ix,:) = evnt(photod_ix,:)-drift_model;
    
    % Correct baseline shift
    for shift_ix = 1:length(bsln_shift_times)
        epoch_idx = floor(bsln_shift_times{shift_ix}(1)*hdr.sample_rate):floor(bsln_shift_times{shift_ix}(2)*hdr.sample_rate);
        epoch_idx(epoch_idx<1) = [];
        evnt(photod_ix,epoch_idx) = evnt(photod_ix,epoch_idx) - bsln_shift_val(shift_ix);
    end
    % zero out drifts
    for zero_ix = 1:length(bsln_times)
        epoch_idx = floor(bsln_times{zero_ix}(1)*hdr.sample_rate):floor(bsln_times{zero_ix}(2)*hdr.sample_rate);
        epoch_idx(epoch_idx<1) = [];
        evnt(photod_ix,epoch_idx) = bsln_val;
    end
    
    % level out stimulus periods
    for stim_ix = 1:length(stim_times)
        epoch_idx = floor(stim_times{stim_ix}(1)*hdr.sample_rate):floor(stim_times{stim_ix}(2)*hdr.sample_rate);
        epoch_idx(epoch_idx<1) = [];
        evnt(photod_ix,epoch_idx) = stim_yval(stim_ix);
    end
    
    % Save corrected data
    out_filename = [SBJ_vars.dirs.preproc SBJ '_evnt_clean',block_suffix,'.mat'];
    save(out_filename, 'evnt', 'hdr', 'ignore_trials', 'photod_ix', 'mic_ix');
    
    %% ========================================================================
    %   Step 5- Parse Event Traces into Behavioral Data
    %  ========================================================================
    if proc.resample_freq~=1000
        error('ERROR!!! SBJ02_behav_parse assumes 1 kHz neural sampling rate!!!\n');
    end
    SU01_behav_parse(SBJ,b_ix,proc.rt_bounds,1,1)
    % Be sure to save the two figures coming from this function!
    %   i.e., SBJ_photodiode_segmentation.fig & SBJ_events.fig
end

%% Plot PSTH
conditions  = 'CNI';
an_opts     = {'PSTH_S_trl2to150_bn20','PSTH_R_trl5to1_bn20'};
plt_opts    = {'ts_trl2to15_errbr_evnt','ts_trl5to1_errbr_evnt'};
plot_ISI    = 1;
fig_vis     = 'on';
save_plots  = 1;
close_plots = 0;
for an_ix = 1%:numel(an_opts)
    SU02_PSTH_ft(SBJ,conditions,proc_id,an_opts{an_ix},plt_opts{an_ix},plot_ISI,fig_vis,save_plots,close_plots);
end

%% ========================================================================
%   Step 6- Create einfo based on ROIs
%  ========================================================================
% look at recon and create spreadsheet of general ROI, WM/GM, etc.
%   save that as tsv
% fn_compile_elec_struct(SBJ,proc_id,'pat')
% fn_compile_elec_struct(SBJ,proc_id,'mni')
% fn_compile_einfo(SBJ,proc_id)

% %% ========================================================================
% %   Step 7- Manually Correct Reaction Times
% %  ========================================================================
% % Send to Rana, get her to create the basic excel sheet
% % double check it!
% 
% % Do this separately for each block
% b_ix = 1;
% if numel(SBJ_vars.raw_file)>1
%     block_suffix = strcat('_',SBJ_vars.block_name{b_ix});
% else
%     block_suffix = SBJ_vars.block_name{b_ix};   % should just be ''
% end
% 
% % Before running this function, open SBJ_events.fig
% % run with default parameters (outlier_thresh=0.1s), not saving plot or trial_info
% outlier_thresh  = 0.15;
% save_plot       = 0;
% save_trial_info = 0;
% SBJ03_RT_manual_adjustments(SBJ,b_ix,outlier_thresh,save_plot,save_trial_info)
% % Check any big discrepancies between auto and manual RTs
% % Fix any RTs manually if necessary
% 
% % Run again after verifying accuracy of manual RTs, saving this time
% save_plot       = 1;
% save_trial_info = 1;
% SBJ03_RT_manual_adjustments(SBJ,b_ix,outlier_thresh,save_plot,save_trial_info)
% 
% %% ========================================================================
% %   Step 8- Preprocess Neural Data
% %  ========================================================================
% SBJ04_preproc(SBJ,proc_id)
% 
% % Second visual cleaning
% load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
% bad_preclean = load(strcat(SBJ_vars.dirs.events,SBJ,'_bob_bad_epochs_preclean',block_suffix,'.mat'));
% bad_at = fn_convert_epochs_full2at(bad_preclean.bad_epochs,SBJ_vars.analysis_time{b_ix},...
%                                             strcat(SBJ_vars.dirs.preproc,SBJ,'_preclean',block_suffix,'.mat'),1);
% load(strcat(root_dir,'PRJ_Error/scripts/utils/cfg_plot.mat'));
% % cfg_plot.artfctdef.visual.artifact = bad_at;
% out = ft_databrowser(cfg_plot,data);
% 
% bad_epochs = out.artfctdef.visual.artifact;
% save(strcat(SBJ_vars.dirs.events,SBJ,'_colin_bad_epochs_preproc.mat'),'-v7.3','bad_epochs');
% 
% %% ========================================================================
% %   Step 9- Reject Bad Trials Based on Behavior and Bob
% %  ========================================================================
% clear data trial_info
% % Load manually corrected trial_info
% ti = {};
% block_lens   = zeros([numel(SBJ_vars.block_name) 1]);
% block_times  = zeros([numel(SBJ_vars.block_name) 1]);
% block_trlcnt = zeros([numel(SBJ_vars.block_name) 1]);
% block_blkcnt = zeros([numel(SBJ_vars.block_name) 1]);
% for b_ix = 1:numel(SBJ_vars.block_name)
%     if numel(SBJ_vars.raw_file)>1
%         block_suffix = strcat('_',SBJ_vars.block_name{b_ix});
%         % Get block length
%         tmp = load(strcat(SBJ_vars.dirs.import,SBJ,'_',...
%             num2str(proc.resample_freq),'hz',block_suffix,'.mat'));
%         block_lens(b_ix) = size(tmp.data.trial{1},2);
%         block_times(b_ix) = tmp.data.time{1}(end);
%     else
%         block_suffix = SBJ_vars.block_name{b_ix};   % should just be ''
%     end
%     ti{b_ix} = load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_manual',block_suffix,'.mat'));
%     block_trlcnt(b_ix) = numel(ti{b_ix}.trial_info.trial_n);
%     block_blkcnt(b_ix) = max(ti{b_ix}.trial_info.block_n);
%     
%     % Add run number
%     ti{b_ix}.trial_info.run_n = ones(size(ti{b_ix}.trial_info.block_n))*b_ix;
% end
% 
% % Combine trial_info structs if necessary
% if numel(ti)<2
%     trial_info = ti{1}.trial_info;
% else
%     trial_info = ti{1}.trial_info;
%     % Add properties of the individual blocks
%     trial_info.run_len        = block_lens;
%     trial_info.run_time       = block_times;
%     trial_info.run_trial_info = ti;     % Keep the original trial_info structs from each run
%     for b_ix = 2:numel(SBJ_vars.block_name)
%         % Concatenate fields that don't need modification
%         %   Not modifying marker_time and onset_time (no idea what those are...)
%         trial_info.word          = vertcat(trial_info.word,ti{b_ix}.trial_info.word);
%         trial_info.color         = vertcat(trial_info.color,ti{b_ix}.trial_info.color);
%         trial_info.trialtype     = vertcat(trial_info.trialtype,ti{b_ix}.trial_info.trialtype);
%         trial_info.blocktype     = vertcat(trial_info.blocktype,ti{b_ix}.trial_info.blocktype);
%         trial_info.response_time = vertcat(trial_info.response_time,ti{b_ix}.trial_info.response_time);
%         trial_info.marker_time   = vertcat(trial_info.marker_time,ti{b_ix}.trial_info.marker_time);
%         trial_info.onset_time    = vertcat(trial_info.onset_time,ti{b_ix}.trial_info.onset_time);
%         trial_info.condition_n   = vertcat(trial_info.condition_n,ti{b_ix}.trial_info.condition_n);
%         trial_info.error         = vertcat(trial_info.error,ti{b_ix}.trial_info.error);
%         trial_info.run_n         = vertcat(trial_info.run_n,ti{b_ix}.trial_info.run_n);
%         
%         % Modify then concatenate counts and indices
%         trial_info.block_n = vertcat(trial_info.block_n,ti{b_ix}.trial_info.block_n+sum(block_blkcnt(1:b_ix-1)));
%         trial_info.trial_n = vertcat(trial_info.trial_n,ti{b_ix}.trial_info.trial_n+sum(block_trlcnt(1:b_ix-1)));
%         trial_info.ignore_trials = horzcat(trial_info.ignore_trials,...
%             ti{b_ix}.trial_info.ignore_trials+sum(block_trlcnt(1:b_ix-1)));
%         
%         trial_info.word_onset = vertcat(trial_info.word_onset,...
%             ti{b_ix}.trial_info.word_onset+sum(block_lens(1:b_ix-1)));
%         trial_info.resp_onset = vertcat(trial_info.resp_onset,...
%             ti{b_ix}.trial_info.resp_onset+sum(block_lens(1:b_ix-1)));
%     end
% end
% clear ti block_lens block_times block_trlcnt block_blkcnt
% 
% % Toss trials based on behavior and cleaning with Bob
% trial_info_clean = SBJ05_reject_behavior(SBJ,trial_info,proc_id);
% 
% %% ========================================================================
% %   Step 10a- Prepare Variance Estimates for Variance-Based Trial Rejection
% %  ========================================================================
% % Load data for visualization
% % load(strcat(preproc_dir,SBJ,'_proc.mat'));
% load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
% 
% % Select channels of interest
% cfg = [];
% cfg.channel = SBJ_vars.ch_lab.ROI;
% data = ft_selectdata(cfg,data);
% 
% % Segment into trials
% if strcmp(proc.event_type,'stim')
%     events = trial_info_clean.word_onset;
% elseif strcmp(proc.event_type,'resp')
%     events = trial_info_clean.resp_onset;
% else
%     error(strcat('ERROR: unknown event_type ',proc.event_type));
% end
% trials = fn_ft_cut_trials_equal_len(data,events,...
%     trial_info_clean.condition_n',proc.trial_lim_s*data.fsample);
% 
% % Compute Derivative
% cfg = [];
% cfg.derivative = 'yes';
% trials_dif = ft_preprocessing(cfg,trials);
% 
% % Compute potential variance limits over trials and channels
% [trial_mat,~] = fn_format_trials_ft2KLA(trials);
% var_mat = std(trial_mat,0,3);
% ch_var_mean = mean(var_mat,2);
% ch_var_thresh = mean(ch_var_mean)+std(ch_var_mean)*proc.var_std_warning_thresh;
% 
% trial_var_mean = mean(var_mat,1);
% trial_var_thresh = mean(trial_var_mean)+std(trial_var_mean)*proc.var_std_warning_thresh;
% 
% [trial_mat_dif,~] = fn_format_trials_ft2KLA(trials_dif);
% var_mat_dif = std(trial_mat_dif,0,3);
% ch_var_mean_dif = mean(var_mat_dif,2);
% ch_var_dif_thresh = mean(ch_var_mean_dif)+std(ch_var_mean_dif)*proc.var_std_warning_thresh;
% 
% trial_var_mean_dif = mean(var_mat_dif,1);
% trial_var_dif_thresh = mean(trial_var_mean_dif)+std(trial_var_mean_dif)*proc.var_std_warning_thresh;
% 
% % Report on potentially bad channels
% bad_var_ch      = trials.label(abs(ch_var_mean) > ch_var_thresh);
% bad_var_dif_ch  = trials.label(abs(ch_var_mean_dif) > ch_var_dif_thresh);
% bad_var_trl     = trial_info_clean.trial_n(abs(trial_var_mean) > trial_var_thresh);
% bad_var_dif_trl = trial_info_clean.trial_n(abs(trial_var_mean_dif) > trial_var_dif_thresh);
% fprintf('==============================================================================================\n');
% fprintf('Simple Variance Rejection:\n');
% fprintf('\tChannel Variance Names: %s\n',bad_var_ch{:});
% fprintf('\tChannel Diff Variance Names: %s\n',bad_var_dif_ch{:});
% fprintf('\tTrial Variance: %i\n',bad_var_trl);
% fprintf('\tTrial Diff Variance: %i\n',bad_var_dif_trl);
% fprintf('==============================================================================================\n');
% 
% % Save results
% var_rej_filename = [SBJ_vars.dirs.events SBJ '_variance_rejection_results.txt'];
% r_file = fopen(var_rej_filename,'a');
% fprintf(r_file,'===============================================================================================\n');
% fprintf(r_file,'Simple Variance Rejection:\n');
% fprintf(r_file,'Run Time: %s\n',datestr(datetime));
% fprintf(r_file,'\tChannel Variance Names: %s\n',bad_var_ch{:});
% fprintf(r_file,'\tChannel Diff Variance Names: %s\n',bad_var_dif_ch{:});
% fprintf(r_file,'\tTrial Variance: %i\n',bad_var_trl);
% fprintf(r_file,'\tTrial Diff Variance: %i\n',bad_var_dif_trl);
% fprintf(r_file,'==============================================================================================\n');
% fclose(r_file);
% 
% %% ========================================================================
% %   Step 10b- Choose Thresholds for Variance-Based Trial Rejection
% %  ========================================================================
% % Visualize data to set limits for variance-based rejection
% cfg_reject = [];
% cfg_reject.method = 'summary';
% ft_rejectvisual(cfg_reject,trials);
% 
% % Visualize Derivative
% ft_rejectvisual(cfg_reject,trials_dif);
% 
% %% ========================================================================
% %   Step 11a- Update Rejection parameters and electrodes based on variance
% %  ========================================================================
% % Comment in evernote note on bad trials and channels!
% % Then the following variables should be written into SBJ_vars:
% 
% % Update SBJ_vars.ch_lab.var_rej field!
% 
% % % Choose thresholds based on plots above
% % artifact_params.std_limit_raw = 7;
% % artifact_params.hard_threshold_raw = 300; % based on maxabs()
% % 
% % artifact_params.std_limit_diff = 7;
% % artifact_params.hard_threshold_diff = 100; % based on maxabs() for trials_dif
% 
% % Re-load SBJ_vars after updating artifact field
% clear SBJ_vars
% clear_cmd = ['clear ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
% eval(clear_cmd); %needed to delete cached version
% eval(SBJ_vars_cmd);
% 
% % Reload data and re-select channels of interest after excluding bad ones
% clear data
% load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
% cfg = [];
% cfg.channel = SBJ_vars.ch_lab.ROI;
% data = ft_selectdata(cfg,data);
% trials = fn_ft_cut_trials_equal_len(data,events,...
%     trial_info_clean.condition_n',proc.trial_lim_s*data.fsample);
% 
% % Compute Derivative
% cfg = [];
% cfg.derivative = 'yes';
% trials_dif = ft_preprocessing(cfg,trials);
% 
% 
% %% ========================================================================
% %   Step 11b- Automatically Reject Bad Trials Based on Variance
% %  ========================================================================
% % Run KLA artifact rejection based on robust variance estimates
% % If too many/few trials are rejected, adjust artifact_params and rerun
% 
% % plot_chans is a cell array of channels to be plotted
% %   {} - empty array indicates skip plotting
% %   'worst' - says plot channels that had > 5 bad trials
% %   full array - list of channel names to plot
% % report [struct] - list of 0/1/2 flags for plotting and reporting options,
% %   0 = no report; 1 = concise report; 2 = verbose report
% %   .hard_thresh: 1=number of trials; 2=print rejected trial_n
% %   .std_thresh: 1=number of trials; 2=print rejected trial_n
% %   .std_plot: 1=plot std distribution with cut off
% plot_ch = {'worst',5};%ft_channelselection({'LPC*','LAC*','RIN*'},data.label);
% report.hard_thresh = 2; % print total rejected and trial numbers
% report.std_thresh  = 1; % print only total rejected
% report.std_plot    = 1; % plot the std distribution and threshold
% 
% trial_info_KLA_clean = SBJ06_reject_artifacts_KLA_report(trials,trial_info_clean,...
%                                         SBJ_vars.artifact_params,plot_ch,report);
% 
% bad_samples = NaN([size(trial_info_KLA_clean.bad_trials.var,1) 2]);
% for t_ix = 1:size(bad_samples,1)
%     bad_samples(t_ix,:) = trials.sampleinfo(trial_info_clean.trial_n==trial_info_KLA_clean.bad_trials.var(t_ix),:);
% end
% cfg = [];
% cfg.continuous = 'no';
% cfg.viewmode = 'vertical';
% cfg.artfctdef.visual.artifact = bad_samples;
% ft_databrowser(cfg, trials);
% 
% ft_databrowser(cfg, trials_dif);
% 
% %% ========================================================================
% %   Step 12- Compile Variance-Based Trial Rejection and Save Results
% %  ========================================================================
% % Re-load SBJ_vars after updating artifact field
% clear SBJ_vars
% eval(clear_cmd); %needed to delete cached version
% eval(SBJ_vars_cmd);
% 
% clear trial_info
% trial_info = trial_info_clean;
% % Document bad trials
% trial_info.bad_trials.variance = SBJ_vars.trial_reject_n';
% trial_info.bad_trials.all = sort([trial_info.bad_trials.all; trial_info.bad_trials.variance]);
% 
% % Remove bad trials
% trial_rejected = ismember(trial_info.trial_n,SBJ_vars.trial_reject_n);
% trial_reject_ix = find(trial_rejected);
% trial_info.block_n(trial_reject_ix) = [];
% trial_info.trial_n(trial_reject_ix) = [];
% trial_info.word(trial_reject_ix) = [];
% trial_info.color(trial_reject_ix) = [];
% trial_info.trialtype(trial_reject_ix) = [];
% trial_info.blocktype(trial_reject_ix) = [];
% trial_info.response_time(trial_reject_ix) = [];
% trial_info.marker_time(trial_reject_ix) = [];
% trial_info.onset_time(trial_reject_ix) = [];
% trial_info.word_onset(trial_reject_ix) = [];
% trial_info.resp_onset(trial_reject_ix) = [];
% trial_info.condition_n(trial_reject_ix) = [];
% trial_info.error(trial_reject_ix) = [];
% 
% % Make sure no response times are in weird float format
% trial_info.resp_onset = round(trial_info.resp_onset);
% 
% save(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
% 
