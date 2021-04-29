%% Pipeline Processing Variables: "main_ft"

% Data Preprocessing
proc_vars.plot_psd      = '1by1';         % type of plot for channel PSDs
proc_vars.resample_yn   = 'yes';
proc_vars.resample_freq = 1000;
proc_vars.demean_yn     = 'yes';
proc_vars.hp_yn         = 'yes';
proc_vars.hp_freq       = 0.5;            % [] skips this step
proc_vars.hp_order      = 4;              % Leaving blank causes instability error, 1 or 2 works
proc_vars.lp_yn         = 'yes';
proc_vars.lp_freq       = 300;            % [] skips this step
proc_vars.notch_type    = 'bandstop';     % method for nothc filtering out line noise

% cleanline   = 'yes';                % Use Tim Mullen's cleanline function
% dft_yn      = 'no';
% bs_yn       = 'no';                % Parameters for this in SBJ_vars

% Behavioral Processing
proc_vars.rt_bounds = [0.5 1.5];          % exclude RTs too early/late 

% Trial Cut Parameteres
proc_vars.event_type    = 'stim';         % 'stim'/'resp': lock trial to these event
proc_vars.trial_lim_s = [-0.25 3];      % data segments (in seconds) to grab around events
proc_vars.RT_std_thresh = 3;              % rejection threshold for RTs

% Varaince-Based Trial Rejection Parameters
proc_vars.var_std_warning_thresh = 3;
