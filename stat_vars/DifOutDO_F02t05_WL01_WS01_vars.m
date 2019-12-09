% Time Parameters
st.evnt_lab = 'F';
st.stat_lim = [0.2 0.5];

% Sliding Window Parameters
st.win_len  = 0.1;
st.win_step = 0.1;

% Model Factors and Levels
st.model_lab = 'DifOutDO';
st.groups    = {'Dif', 'Out'};
st.anova_term = 'interaction';
st.alpha     = 0.05;
st.n_boots   = 1000;

% ANOVA Parameters
st.regress_rt = 0;    % regrees reaction time off before running ANOVA

% RT Correlation Parameters
st.rt_corr = 0;
cfg_rt = [];
cfg_rt.parameter        = 'powspctrm';
cfg_rt.statistic        = 'ft_statfun_correlationT';
cfg_rt.method           = 'montecarlo';
cfg_rt.numrandomization = st.n_boots;
cfg_rt.correctm         = 'cluster';
cfg_rt.clusteralpha     = 0.05;   %threshold for a single comparison (time point) to be included in the clust
cfg_rt.clusterstatistic = 'maxsum';
cfg_rt.clustertail      = 0;
cfg_rt.tail             = 0; %two sided
cfg_rt.correcttail      = 'alpha'; %correct the .alpha for two-tailed test (/2)
cfg_rt.computestat      = 'yes';
cfg_rt.computeprob      = 'yes';
cfg_rt.alpha            = 0.05;
cfg_rt.neighbours       = [];

