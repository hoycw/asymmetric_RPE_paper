event_type  = 'fb';           % event around which to cut trials
% trial_lim_s will NOT be full of data! the first and last t_ftimwin/2 epochs will be NaNs
trial_lim_s = [-0.25 1.21];      % window in SEC for cutting trials
%plt_lim     = [-0.2 1.5];         % window to plot this data
demean_yn   = 'no';             % z-score for HFA instead
bsln_evnt   = 'stim';
bsln_type   = 'zboot';
bsln_lim    = [-0.25 -0.05];    % window in SEC for baseline correction

% HFA Calculations
HFA_type = 'multiband';                 % b = broadband = 70-150 Hz
foi_center  = [70:10:150];
octave      = 3/4;              % Frequency resolution
foi_min     = 2^(-octave/2)*foi_center;
foi_max     = 2^(octave/2)*foi_center;
foi         = (foi_min+foi_max)/2;
delta_freq  = foi_max-foi_min;
delta_time  = 0.1;
n_taper_all = max(1,round(delta_freq.*delta_time-1));   %number of tapers for each frequency
foi_center  = round(foi_center*10)/10;          %convert to float?
delta_freq_true = (n_taper_all+1)./delta_time; % total bandwidth around

cfg_hfa = [];
cfg_hfa.output       = 'pow';
cfg_hfa.channel      = 'all';
cfg_hfa.method       = 'mtmconvol';
cfg_hfa.taper        = 'dpss';
cfg_hfa.tapsmofrq    = delta_freq_true./2;                  %ft wants half bandwidth around the foi
cfg_hfa.keeptapers   = 'no';
cfg_hfa.pad          = 'maxperlen';                         %add time on either side of window
cfg_hfa.padtype      = 'zero';
cfg_hfa.foi          = foi_center;                          % analysis 2 to 30 Hz in steps of 2 Hz 
cfg_hfa.t_ftimwin    = ones(length(cfg_hfa.foi),1).*delta_time;    % length of time window; 0.5 sec, could be n_cycles./foi for n_cylces per win
cfg_hfa.toi          = 'all';%-buff_lim(1):0.1:1.5;         % time window centers
cfg_hfa.keeptrials   = 'yes';                               % must be 'yes' for stats
% cfg.t_ftimwin    = ones(1,length(cfg.tapsmofrq))*delta_time;

% Log Tranform
log_yn = 0;

% Outlier Rejection
outlier_std_lim = 6;

% Cleaning up power time series for plotting
smooth_pow_ts = 0;
lp_yn       = 'no';
lp_freq     = 10;
hp_yn       = 'no';
hp_freq     = 0.5;

% Resampling
resample_ts   = 0;
% resample_freq = 250;

% Stats parameters
stat_lim    = [0 1];            % window in SEC for stats
n_boots     = 1000;             % Repetitions for non-parametric stats

cfg_stat = [];
cfg_stat.latency          = stat_lim;
cfg_stat.channel          = 'all';
cfg_stat.parameter        = 'powspctrm';
cfg_stat.method           = 'montecarlo';
cfg_stat.statistic        = 'ft_statfun_indepsamplesT';
cfg_stat.correctm         = 'cluster';
cfg_stat.clusteralpha     = 0.05;   %threshold for a single comparison (time point) to be included in the clust
cfg_stat.clusterstatistic = 'maxsum';
cfg_stat.clustertail      = 0;
cfg_stat.tail             = 0; %two sided
cfg_stat.correcttail      = 'alpha'; %correct the .alpha for two-tailed test (/2)
cfg_stat.alpha            = 0.05;
cfg_stat.numrandomization = n_boots;
cfg_stat.neighbours       = [];%neighbors;
% cfg_stat.minnbchan        = 0;
cfg_stat.ivar             = 1;  %row of design matrix containing independent variable
% cfg_stat.uvar             = 2;  %row containing dependent variable, not needed for indepsamp

