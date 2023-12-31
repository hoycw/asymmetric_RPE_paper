% Stat Parameters
st.an_style    = 'lme';
st.stat_cond   = 'DifFB';
st.stat_lim    = [-0.05 0.05];            % window in SEC for stats
st.measure     = 'erp_mean';             % {'ts', 'p2p', 'mean', 'erp_mean'}
st.n_boots     = 1000;             % Repetitions for non-parametric stats
st.alpha       = 0.05;
st.mcp_method  = 'FDR';

% Peak time selection parameters
st.pk_an_id    = 'ERP_Pz_F2t1_dm2t0_fl05t20';            % an_id from which to get peak time
st.pk_cond_grp = 'All';                         % set of conditions in ERP analysis to load
st.pk_erp_cond = 'All';                         % set of conditions to select within ERP analysis to get peak info
st.pk_lim      = [0.25 0.5];                    % Window to search for peak
st.pk_sign     = 1;                            % Sign of peak to find

