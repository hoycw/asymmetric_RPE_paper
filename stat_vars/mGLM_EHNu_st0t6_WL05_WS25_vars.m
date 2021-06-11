% Stat Parameters
st.an_style    = 'mGLM';
st.stat_cond   = 'EHNu';
st.stat_lim    = [0 0.6];            % window in SEC for stats

% Sliding Window Parameters
st.win_len  = 0.05;
st.win_step = 0.025;
st.win_center = [st.stat_lim(1)+st.win_len/2:st.win_step:st.stat_lim(2)-st.win_len/2];

st.measure     = 'ts';             % {'ts', 'p2p', 'mean'}
st.n_boots     = 1000;             % Repetitions for non-parametric stats
st.alpha       = 0.05;
st.mcp_method  = 'FDR';
st.test_tail   = 2;                 % one-sided (1) or two-sided (2) tests

% RT correlation and regression
st.rt_corr     = 0;
st.regress_rt  = 0;

if st.rt_corr
    cfg_rt = [];
    cfg_rt.parameter        = 'powspctrm';
    cfg_rt.statistic        = 'ft_statfun_correlationT';
    cfg_rt.method           = 'montecarlo';
    cfg_rt.numrandomization = st.n_boots;
    cfg_rt.alpha            = st.alpha;
%     cfg_rt.correctm         = 'cluster';
%     cfg_rt.clusteralpha     = 0.05;   %threshold for a single comparison (time point) to be included in the clust
%     cfg_rt.clusterstatistic = 'maxsum';
%     cfg_rt.clustertail      = 0;
%     cfg_rt.tail             = 0; %two sided
%     cfg_rt.correcttail      = 'alpha'; %correct the .alpha for two-tailed test (/2)
%     cfg_rt.computestat      = 'yes';
%     cfg_rt.computeprob      = 'yes';
%     cfg_rt.neighbours       = [];
end
