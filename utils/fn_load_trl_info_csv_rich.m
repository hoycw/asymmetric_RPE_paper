function [trl_info] = fn_load_trl_info_csv_rich(csv_filename,ignore_trials)
error('this was for the rich version with unnecessary variables');
%% Load trial info csv exported from python
%   Also converts from python 0-based to MATLAB 1-based indexing
% INPUTS:
%   csv_filename [str] - subject ID to load
%   ignore_trials [int array] - list of trial numbers to be removed (MATLAB indexing)
% OUTPUTS:
%   trl_info [struct] - contains many fields of info on the trials

% Input format:
%   Total_Trial [int] - trial # from 1:n_trials*n_blocks (+training+examples)
%   Block [int] - block #
%   Condition [str] - 'easy' or 'hard'
%   Hit [int] - 0/1 if this trial was correct
%   RT [float] - reaction time on this trial
%   Reversal [int] - 0/1 if this trial outcome was different from last
%       (reversed staircase direction)
%   Timestamp [float] - time of outcome logging (since start of experiment)
%   Tolerance [float] - tolerance around target interval for this trial
%   Trial [int] - trial # within a block (e.g., 1:n/trials/block_n)
%   Score [int] - score obtained on this trial (+/- 100)
%   ITI [float] - time since timestamp n-1, will be 0 for first trial in block, training data
%   ITI Type [str] - ['short','medium', or 'long']
%       only have medium in cases with 3 (not 4) ITIs (e.g., IR57)
%   PE [str] - True/False logical for last trial being an error
%   Score Total [int] -  total score accrued up to this trial
%   Hit Total [int] - total number of hits accrued so far
%   Accuracy [float] - accuracy on this block (entire block!)

fields = {'trial_n','block','cond','hit','rt','reversal','time',...
            'tolerance','block_trial','score','ITI','ITI_type','post_err',...
            'score_total','hit_total','block_acc'};
formatspec = '%d%d%s%d%f%d%f%f%d%d%f%s%s%d%d%f';
csv_file = fopen(csv_filename);
csv = textscan(csv_file,formatspec,'Delimiter',',','HeaderLines',1);
fclose(csv_file);

for ix = 1:numel(fields)
    switch fields{ix}
        case 'Total_Trial'
            trl_info.trial_n = str2double(csv.textdata(2:end,ix))+1;
        case 'Block'
            trl_info.block = str2double(csv.textdata(2:end,ix))+1;
        case 'Condition'
            trl_info.cond = csv.textdata(2:end,ix);
        case 'Hit'
            trl_info.hit = str2double(csv.textdata(2:end,ix));
        case 'RT'
            trl_info.rt = str2double(csv.textdata(2:end,ix));
        case 'Reversal'
            trl_info.reversal = str2double(csv.textdata(2:end,ix));
        case 'Timestamp'
            trl_info.time = str2double(csv.textdata(2:end,ix));
        case 'Tolerance'
            trl_info.tolerance = str2double(csv.textdata(2:end,ix));
        case 'Trial'
            trl_info.block_trial = str2double(csv.textdata(2:end,ix))+1;
        case 'Score'
            trl_info.score = str2double(csv.textdata(2:end,ix));
        case 'Score Total'
            trl_info.score_total = str2double(csv.textdata(2:end,ix));
        case 'Hit Total'
            trl_info.hit_total = str2double(csv.textdata(2:end,ix));
        case 'ITI'
            trl_info.ITI = str2double(csv.textdata(2:end,ix));
        case 'ITI type'
            trl_info.ITI_type = csv.textdata(2:end,ix);
        case 'PE'
            trl_info.post_err = false(size(csv.textdata,1)-1,1);
            for t_ix = 2:size(csv.textdata,1)
                trl_info.post_err(t_ix-1) = strcmpi(csv.textdata(t_ix,ix), 'True');
            end
        case 'Accuracy'
            trl_info.block_acc = str2double(csv.textdata(2:end,ix));
        otherwise
            error(['Unknown field encountered: ' fields{ix}]);
    end
end

% Remove ignored trials
if ~isempty(ignore_trials)
    new_fields = fieldnames(trl_info);
    n_trials = numel(trl_info.trial_n);
    for f_ix = 1:numel(new_fields)
        if eval(['numel(trl_info.' new_fields{f_ix} ')==n_trials'])
            eval(['trl_info.' new_fields{f_ix} '(ignore_trials) = [];']);
        end
    end
end
