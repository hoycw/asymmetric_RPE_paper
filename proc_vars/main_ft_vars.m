%% Pipeline Processing Variables: "main_ft"

% Data Preprocessing
proc.plot_psd      = '1by1';         % type of plot for channel PSDs
proc.resample_yn   = 'yes';
proc.resample_freq = 1000;
proc.demean_yn     = 'yes';
proc.hp_yn         = 'yes';
proc.hp_freq       = 0.5;            % [] skips this step
proc.hp_order      = 4;              % Leaving blank causes instability error, 1 or 2 works
proc.lp_yn         = 'yes';
proc.lp_freq       = 300;            % [] skips this step
proc.notch_type    = 'bandstop';     % method for nothc filtering out line noise

% cleanline   = 'yes';                % Use Tim Mullen's cleanline function
% dft_yn      = 'no';
% bs_yn       = 'no';                % Parameters for this in SBJ

% Behavioral Processing
proc.rt_bounds = [0.5 1.5];          % exclude RTs too early/late 

% Trial Cut Parameteres
proc.event_type    = 'stim';         % 'stim'/'resp': lock trial to these event
proc.trial_lim_s = [-0.25 3];      % data segments (in seconds) to grab around events
proc.RT_std_thresh = 3;              % rejection threshold for RTs

% Varaince-Based Trial Rejection Parameters
proc.var_std_warning_thresh = 3;
