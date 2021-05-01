function [combined] = fn_combine_trl_info(ti_cell)
%% Combine trl_info structs from multiple runs into a single struct
% INPUTS:
%   trl_info [cell array] - should have multiple trl_info structs
% OUTPUTS:
%   combined [struct] - single trl_info struct with combined data

ti = {};
block_lens   = zeros([numel(SBJ_vars.block_name) 1]);
block_times  = zeros([numel(SBJ_vars.block_name) 1]);
block_trlcnt = zeros([numel(SBJ_vars.block_name) 1]);
block_blkcnt = zeros([numel(SBJ_vars.block_name) 1]);
for b_ix = 1:numel(SBJ_vars.block_name)
    if numel(SBJ_vars.raw_file)>1
        block_suffix = strcat('_',SBJ_vars.block_name{b_ix});
        % Get block length
        tmp = load(strcat(SBJ_vars.dirs.import,SBJ,'_',...
            num2str(proc.resample_freq),'hz',block_suffix,'.mat'));
        block_lens(b_ix) = size(tmp.data.trial{1},2);
        block_times(b_ix) = tmp.data.time{1}(end);
    else
        block_suffix = SBJ_vars.block_name{b_ix};   % should just be ''
    end
    ti{b_ix} = load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_manual',block_suffix,'.mat'));
    block_trlcnt(b_ix) = numel(ti{b_ix}.trial_info.trial_n);
    block_blkcnt(b_ix) = max(ti{b_ix}.trial_info.block_n);
    
    % Add run number
    ti{b_ix}.trial_info.run_n = ones(size(ti{b_ix}.trial_info.block_n))*b_ix;
end

combined = trl_info{1};
% Add properties of the individual blocks
trl_info.run_len        = block_lens;
trl_info.run_time       = block_times;
trl_info.run_trl_info   = ti_cell;     % Keep the original trl_info structs from each run
for b_ix = 2:numel(SBJ_vars.block_name)
    % Concatenate fields that don't need modification
    %   Not modifying marker_time and onset_time (no idea what those are...)
    trl_info.word          = vertcat(trl_info.word,ti_cell{b_ix}.trl_info.word);
    trl_info.color         = vertcat(trl_info.color,ti_cell{b_ix}.trl_info.color);
    trl_info.trialtype     = vertcat(trl_info.trialtype,ti_cell{b_ix}.trl_info.trialtype);
    trl_info.blocktype     = vertcat(trl_info.blocktype,ti_cell{b_ix}.trl_info.blocktype);
    trl_info.response_time = vertcat(trl_info.response_time,ti_cell{b_ix}.trl_info.response_time);
    trl_info.marker_time   = vertcat(trl_info.marker_time,ti_cell{b_ix}.trl_info.marker_time);
    trl_info.onset_time    = vertcat(trl_info.onset_time,ti_cell{b_ix}.trl_info.onset_time);
    trl_info.error         = vertcat(trl_info.error,ti_cell{b_ix}.trl_info.error);
    trl_info.run_n         = vertcat(trl_info.run_n,ti_cell{b_ix}.trl_info.run_n);
    trl_info.condition_n   = horzcat(trl_info.condition_n,ti_cell{b_ix}.trl_info.condition_n);
    
    % Modify then concatenate counts and indices
    trl_info.block_n = vertcat(trl_info.block_n,ti_cell{b_ix}.trl_info.block_n+sum(block_blkcnt(1:b_ix-1)));
    trl_info.trial_n = vertcat(trl_info.trial_n,ti_cell{b_ix}.trl_info.trial_n+sum(block_trlcnt(1:b_ix-1)));
    trl_info.ignore_trials = horzcat(trl_info.ignore_trials,...
        ti_cell{b_ix}.trl_info.ignore_trials+sum(block_trlcnt(1:b_ix-1)));
    
    trl_info.word_onset = vertcat(trl_info.word_onset,...
        ti_cell{b_ix}.trl_info.word_onset+sum(block_lens(1:b_ix-1)));
    trl_info.resp_onset = vertcat(trl_info.resp_onset,...
        ti_cell{b_ix}.trl_info.resp_onset+sum(block_lens(1:b_ix-1)));
end


end