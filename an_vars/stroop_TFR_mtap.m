event_type  = 'stim';           % event around which to cut trials
% trial_lim_s will be expanded in SBJ09a by t_ftimwin/2 on front and back to avoid NaNs within real trial_lim_s
trial_lim_s = [-0.5 2];         % window in SEC for cutting trials
%plt_lim     = [-0.5 2];         % window to plot this data
demean_yn   = 'no';
bsln_evnt   = 'stim';
bsln_type   = 'relchange';
bsln_lim    = [-0.3 -0.1];    % window in SEC for baseline correction

% TFR Calculations
foi_center  = [2:0.5:6 7:16 18 20 25 30 40 50 60 80 100 125 150 175 200 250];%2.^[2:1/4:8]; % Center frequencies
octave      = 3/4;              % Frequency resolution
foi_min     = 2^(-octave/2)*foi_center;
foi_max     = 2^(octave/2)*foi_center;
foi         = (foi_min+foi_max)/2;
delta_freq  = foi_max-foi_min;
delta_time  = 0.5;
n_taper_all = max(1,round(delta_freq.*delta_time-1));   %number of tapers for each frequency
foi_center  = round(foi_center*10)/10;          %convert to float?
delta_freq_true = (n_taper_all+1)./delta_time; % total bandwidth around

cfg_tfr = [];
cfg_tfr.output       = 'pow';
cfg_tfr.channel      = 'all';
cfg_tfr.method       = 'mtmconvol';
cfg_tfr.taper        = 'dpss';
cfg_tfr.tapsmofrq    = delta_freq_true./2;                  %ft wants half bandwidth around the foi
cfg_tfr.keeptapers   = 'no';
cfg_tfr.pad          = 'maxperlen';                         %add time on either side of window
cfg_tfr.padtype      = 'zero';
cfg_tfr.foi          = foi_center;%2:5:150;                 % analysis 2 to 30 Hz in steps of 2 Hz 
cfg_tfr.t_ftimwin    = ones(length(cfg_tfr.foi),1).*0.5;    % length of time window; 0.5 sec, could be n_cycles./foi for n_cylces per win
cfg_tfr.toi          = 'all';%-buff_lim(1):0.1:1.5;         % time window centers
cfg_tfr.keeptrials   = 'yes';                               % must be 'yes' for stats
% cfg.t_ftimwin    = ones(1,length(cfg.tapsmofrq))*delta_time;

% Stats parameters
stat_lim    = [0 2];            % window in SEC for stats
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

