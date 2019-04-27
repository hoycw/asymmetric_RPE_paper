event_type  = 'stim';           % event around which to cut trials
% trial_lim_s will NOT be full of data! the first and last t_ftimwin/2 epochs will be NaNs
trial_lim_s = [-0.25 1.51];      % window in SEC for cutting trials
demean_yn   = 'no';             % z-score for HFA instead
bsln_evnt   = 'stim';
bsln_type   = 'zboot';
bsln_lim    = [-0.25 -0.05];    % window in SEC for baseline correction

% HFA Calculations
HFA_type   = 'hilbert';
foi_lim = [70 150]; % min and max of desired frequencies
n_foi   = 8;
min_exp = log(foi_lim(1))/log(2); % that returns the exponents
max_exp = log(foi_lim(2))/log(2);
fois    = 2.^[linspace(min_exp,max_exp,n_foi)];
foi_bws = fn_semilog_bws(fois);     % semilog bandwidth spacing to match Erik Edwards & Chang lab
bp_lim  = zeros([numel(fois) 2]);
for f = 1:numel(fois)
    bp_lim(f,:) = fn_freq_lim_from_CFBW(fois(f), foi_bws(f));
end

cfg_hfa = [];
cfg_hfa.hilbert  = 'abs';
cfg_hfa.bpfilter = 'yes';
cfg_hfa.bpfreq   = [];      % to be filled by looping through foi_center
cfg_hfa.channel  = 'all';

% Outlier Rejection
% outlier_std_lim = 6;

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
stat_lim    = [0 1.5];            % window in SEC for stats
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

