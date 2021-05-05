function [bhv_trim] = fn_select_bhv(bhv, cond_idx)
%% Select trials within bhv struct (all non-zero entries in cond_idx)

bhv_trim = bhv;
bhv_fields = fieldnames(bhv);
orig_n_trials = numel(bhv.trl_n);
for f_ix = 1:numel(bhv_fields)
    if numel(bhv.(bhv_fields{f_ix}))==orig_n_trials
        bhv_trim.(bhv_fields{f_ix}) = bhv.(bhv_fields{f_ix})(cond_idx~=0);
    end
end

end