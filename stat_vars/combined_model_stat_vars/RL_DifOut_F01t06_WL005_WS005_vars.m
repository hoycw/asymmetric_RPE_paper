% Time Parameters
st.evnt_lab = 'F';
st.stat_lim = [0.1 0.6];

% Sliding Window Parameters
st.win_len  = 0.05;
st.win_step = 0.05;

% Model Factors and Levels
st.model_lab  = 'pWinPEus';
st.regressors = {'pWin', 'sPE', 'uPE', 'offset'};
st.add_offset = 1;
st.trial_cond = {'DifOut'};% needs 'DifFB'
st.alpha      = 0.05;
st.n_boots    = 1000;
%st.mcp_method  = 'FDR';

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

