%% Preprocessing Pipeline
% This script should be run in sections. Functions/scripts with the SBJ##
% prefix can be run automatically, and all other sections should be
% manually editted for each dataset.

% Set Up Directories
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath(genpath([root_dir 'PRJ_Error/scripts/']));
addpath(genpath([root_dir 'PRJ_Error/scripts/utils/']));
addpath(ft_dir);
ft_defaults

%% Step 0 - Processing Variables
% SBJ = 'IR';
proc_id = 'main_ft';
eval(['run ' root_dir 'PRJ_Error/scripts/proc_vars/' proc_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m']);

%% ======================================================================== 
%   Step 1- Quick Import and Processing for Data Cleaning/Inspection
%  ========================================================================
% FILE TOO BIG, RUNNING THIS VIA SGE
% SBJ00a_cleaning_prep(SBJ,proc.plot_psd);

%% ========================================================================
%   Step 2- Initial Viewing of Raw Data (colored by QA)
%  ========================================================================
block_ix    = 1;
keep_db_out = 1;
reorder     = {};   % {} for alphabetical
browser_out = SBJ00b_view_preclean(SBJ,block_ix,keep_db_out,'reorder',reorder);

% Save out the bad_epochs from the preprocessed data
bad_epochs = browser_out.artfctdef.visual.artifact;
tiny_bad = find(diff(bad_epochs,1,2)<10);
if ~isempty(tiny_bad)
    warning([num2str(numel(tiny_bad)) ' tiny bad epochs detected:\n']);
    disp(bad_epochs(tiny_bad,:));
    bad_epochs(tiny_bad,:) = [];
end
save(strcat(SBJ_vars.dirs.events,SBJ,'_bad_epochs_preclean.mat'),'-v7.3','bad_epochs');

%% ========================================================================
%   Step 3- Import Data, Resample, and Save Individual Data Types
%  ========================================================================
%   Run this after rejecting bad channels from preclean viewing
%   If file is too big, run this via SGE

SBJ01a_import_data(SBJ,proc_id);

nlx_align_save_it = 1;
if isfield(SBJ_vars.dirs,'nlx')
    SBJ01b_align_nlx_evnt(SBJ,proc_id,block_ix,nlx_align_save_it);
end

%% ========================================================================
%   Step 3- Preprocess Neural Data
%  ========================================================================
SBJ02_preproc(SBJ,proc_id);

%% ========================================================================
%   Step 4- Second visual cleaning
%  ========================================================================
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
% Load bad_epochs from preclean data and adjust to analysis_time
preclean_ep_at = fn_compile_epochs_full2at(SBJ,proc_id);

% Plot data with bad_epochs highlighted
load(strcat(root_dir,'PRJ_Error/scripts/utils/cfg_plot.mat'));
% If you want to see preclean bad_epochs:
cfg_plot.artfctdef.visual.artifact = preclean_ep_at;
if isfield(data,'sampleinfo')
    data = rmfield(data,'sampleinfo');
end
out = ft_databrowser(cfg_plot,data);

% Save out the bad_epochs from the preprocessed data
bad_epochs = out.artfctdef.visual.artifact;
tiny_bad = find(diff(bad_epochs,1,2)<10);
if ~isempty(tiny_bad)
    warning('Tiny bad epochs detected:\n');
    disp(bad_epochs(tiny_bad,:));
    bad_epochs(tiny_bad,:) = [];
end
save(strcat(SBJ_vars.dirs.events,SBJ,'_bad_epochs_preproc.mat'),'-v7.3','bad_epochs');

%% ========================================================================
%   Step 5a- Manually Clean Photodiode Trace: Load & Plot
%  ========================================================================
% Load data
trl_info     = cell(size(SBJ_vars.block_name));
trl_info_cln = cell(size(SBJ_vars.block_name));
for b_ix = 1:numel(SBJ_vars.block_name)
    % Create a block suffix in cases with more than one recording block
    if numel(SBJ_vars.raw_file)==1 || isfield(SBJ_vars.dirs,'nlx')
        block_suffix = SBJ_vars.block_name{b_ix};   % should just be ''
    else
        block_suffix = strcat('_',SBJ_vars.block_name{b_ix});
    end
    if SBJ_vars.low_srate(b_ix)~=0
        evnt_srate = SBJ_vars.low_srate(b_ix);
    else
        evnt_srate = proc.resample_freq;
    end
    evnt_fname = strcat(SBJ_vars.dirs.import,SBJ,'_evnt_',num2str(evnt_srate),'hz',block_suffix,'.mat');
    load(evnt_fname);
    
    % Plot event channels
    plot(evnt.time{1}, evnt.trial{1});
    
%     % Save .edf for pd-parser
%     evnt_hdr  = ft_fetch_header(evnt);
%     evnt_edf_fname = strcat(SBJ_vars.dirs.import,SBJ,'_evnt_',num2str(evnt_srate),'hz',block_suffix,'.edf');
%     ft_write_data(evnt_edf_fname,evnt.trial{1},'header',evnt_hdr);
    
    %% ========================================================================
    %   Step 5b- Manually Clean Photodiode Trace: Mark Sections to Correct
    %  ========================================================================
    % Create correction times and values in a separate file in ~/PRJ_Stroop/scripts/SBJ_evnt_clean/
    SBJ_evnt_clean_cmd = ['run ' root_dir 'PRJ_Error/scripts/SBJ_evnt_clean/' SBJ '_evnt_clean_params',block_suffix,'.m'];
    eval(SBJ_evnt_clean_cmd);
    
    %% ========================================================================
    %   Step 5c- Manually Clean Photodiode Trace: Apply Corrections
    %  ========================================================================
    % Correct baseline shift
    for shift_ix = 1:length(bsln_shift_times)
        epoch_idx = floor(bsln_shift_times{shift_ix}(1)*evnt.fsample):floor(bsln_shift_times{shift_ix}(2)*evnt.fsample);
        epoch_idx(epoch_idx<1) = [];
        evnt.trial{1}(epoch_idx) = evnt.trial{1}(epoch_idx) - bsln_shift_val(shift_ix);
    end
    % zero out drifts
    for zero_ix = 1:length(bsln_times)
        epoch_idx = floor(bsln_times{zero_ix}(1)*evnt.fsample):floor(bsln_times{zero_ix}(2)*evnt.fsample);
        epoch_idx(epoch_idx<1) = [];
        evnt.trial{1}(epoch_idx) = bsln_val;
    end
    
    % level out stimulus periods
    for stim_ix = 1:length(stim_times)
        epoch_idx = floor(stim_times{stim_ix}(1)*evnt.fsample):floor(stim_times{stim_ix}(2)*evnt.fsample);
        epoch_idx(epoch_idx<1) = [];
        evnt.trial{1}(epoch_idx) = stim_yval(stim_ix);
    end
    
    %% Save corrected data
    out_fname = [SBJ_vars.dirs.preproc SBJ '_evnt_clean',block_suffix,'.mat'];
    save(out_fname, 'evnt', 'ignore_trials');
    
    %% ========================================================================
    %   Step 6- Parse Event Traces into Behavioral Data
    %  ========================================================================
    [~] = SBJ03_behav_parse(SBJ,b_ix,proc_id,1,1);
end

%% ========================================================================
%   Step 7- Reject Bad Trials (Behavior, Visual Cleaning) and Compile Runs
%  ========================================================================
clear trl_info
trl_info = SBJ04_compile_clean_behavior(SBJ,proc_id,1);

%% ========================================================================
%   Step 8- Create elec files from recon
%  ========================================================================
% fn_compile_elec_struct(SBJ,'main_ft','pat','',1);
% fn_compile_elec_struct(SBJ,'main_ft','mni','v',1);
% fn_save_elec_atlas(SBJ,'main_ft','pat','','DK',1);
% fn_save_elec_atlas(SBJ,'main_ft','pat','','Dx',1);
% tissue compartments...

%% %% ========================================================================
% %   Step 9a- Prepare Variance Estimates for Variance-Based Trial Rejection
% %  ========================================================================
% % Load data for visualization
% clear data
% load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
% 
% % Select channels of interest
% cfg = [];
% cfg.channel = SBJ_vars.ch_lab.ROI;
% data = ft_selectdata(cfg,data);
% 
% % Segment into trials
% if strcmp(proc.event_type,'stim')
%     events = trl_info_cln.trl_onset;
% elseif strcmp(proc.event_type,'resp')
%     events = trl_info_cln.rsp_onset;
% else
%     error(stract('ERROR: unknown event_type ',proc.event_type));
% end
% trials = fn_ft_cut_trials_equal_len(data,events,...
%     trl_info_cln.condition_n',proc.trial_lim_s*data.fsample);
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
% bad_var_trl     = trl_info_cln.trl_n(abs(trial_var_mean) > trial_var_thresh);
% bad_var_dif_trl = trl_info_cln.trl_n(abs(trial_var_mean_dif) > trial_var_dif_thresh);
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
% fprintf(r_file,'==========================================================h====================================\n');
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
% %   Step 9b- Choose Thresholds for Variance-Based Trial Rejection
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
% %   Step 10a- Update Rejection parameters and electrodes based on variance
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
% clear_cmd = ['clear /home/knight/hoycw/PRJ_Error/scripts/SBJ_vars/' SBJ '_vars.m'];
% eval(clear_cmd); %needed to delete cached version
% eval(SBJ_vars_cmd);
% 
% % Reload data and re-select channels of interest after excluding bad ones
% clear data
% load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
% cfg = [];
% cfg.channel = SBJ_vars.ch_lab.ROI;
% data = ft_selectdata(cfg,data);
% 
% trials = fn_ft_cut_trials_equal_len(data,events,...
%     trl_info_cln.condition_n',proc.trial_lim_s*data.fsample);
% 
% 
% %% ========================================================================
% %   Step 10b- Automatically Reject Bad Trials Based on Variance
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
% %   .hard_thresh: 1=number of trials; 2=print rejected trl_n
% %   .std_thresh: 1=number of trials; 2=print rejected trl_n
% %   .std_plot: 1=plot std distribution with cut off
% plot_ch = {'worst',5};%ft_channelselection({'LPC*','LAC*','RIN*'},data.label);
% report.hard_thresh = 2; % print total rejected and trial numbers
% report.std_thresh  = 1; % print only total rejected
% report.std_plot    = 1; % plot the std distribution and threshold
% 
% trl_info_KLA_clean = SBJ06_reject_artifacts_KLA_report(trials,trl_info_cln,...
%                                         SBJ_vars.artifact_params,plot_ch,report);
% 
% bad_samples = NaN([size(trl_info_KLA_clean.bad_trials.var,1) 2]);
% for t_ix = 1:size(bad_samples,1)
%     bad_samples(t_ix,:) = trials.sampleinfo(trl_info_cln.trl_n==trl_info_KLA_clean.bad_trials.var(t_ix),:);
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
% %   Step 11- Compile Variance-Based Trial Rejection and Save Results
% %  ========================================================================
% % Re-load SBJ_vars after updating artifact field
% clear SBJ_vars
% eval(clear_cmd); %needed to delete cached version
% eval(SBJ_vars_cmd);
% 
% clear trl_info
% trl_info = trl_info_cln;
% % Document bad trials
% trl_info.bad_trials.variance = SBJ_vars.trial_reject_n';
% trl_info.bad_trials.all = sort([trl_info.bad_trials.all; trl_info.bad_trials.variance]);
% 
% % Remove bad trials
% trial_rejected = ismember(trl_info.trl_n,SBJ_vars.trial_reject_n);
% trial_reject_ix = find(trial_rejected);
% trl_info.block_n(trial_reject_ix) = [];
% trl_info.trl_n(trial_reject_ix) = [];
% trl_info.word(trial_reject_ix) = [];
% trl_info.color(trial_reject_ix) = [];
% trl_info.trialtype(trial_reject_ix) = [];
% trl_info.blocktype(trial_reject_ix) = [];
% trl_info.response_time(trial_reject_ix) = [];
% trl_info.marker_time(trial_reject_ix) = [];
% trl_info.onset_time(trial_reject_ix) = [];
% trl_info.trl_onset(trial_reject_ix) = [];
% trl_info.rsp_onset(trial_reject_ix) = [];
% trl_info.condition_n(trial_reject_ix) = [];
% trl_info.error(trial_reject_ix) = [];


