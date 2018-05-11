function [trl_info] = fn_load_trl_info_csv(csv_filename,ignore_trials)
%% Load trial info csv exported from python
%   Also converts from python 0-based to MATLAB 1-based indexing
% INPUTS:
%   csv_filename [str] - subject ID to load
%   ignore_trials [int array] - list of trial numbers to be removed (MATLAB indexing)
% OUTPUTS:
%   trl_info [struct] - contains many fields of info on the trials

% CSV format:
%   Total_Trial [int] - trial # from 1:n_trials*n_blocks (+training+examples)
%   Block [int] - block #
%   Condition [str] - 'easy' or 'hard'
%   Hit [int] - 0/1 if this trial was correct
%   RT [float] - reaction time on this trial
%   Timestamp [float] - time of outcome logging (since start of experiment)
%   Tolerance [float] - tolerance around target interval for this trial
%   Trial [int] - trial # within a block (e.g., 1:n/trials/block_n)
%   Score [int] - score obtained on this trial (+/- 100)
%   ITI [float] - time since timestamp n-1, will be 0 for first trial in block, training data
%   ITI Type [str] - ['short','medium', or 'long']
%       only have medium in cases with 3 (not 4) ITIs (e.g., IR57)

fields = {'trl_n','blk','cond','hit','rt','time',...
            'tol','blk_trl','score','ITI','ITI_type'};
formatspec = '%d%d%s%d%f%f%f%d%d%f%s';
csv_file = fopen(csv_filename);
csv = textscan(csv_file,formatspec,'Delimiter',',','HeaderLines',1);
fclose(csv_file);

for ix = 1:numel(fields)
    trl_info.(fields{ix}) = csv{ix};
end

% Convert from python to MATLAB indexing
trl_info.trl_n   = trl_info.trl_n+1;
trl_info.blk     = trl_info.blk+1;
trl_info.blk_trl = trl_info.blk_trl+1;

% Remove ignored trials
if ~isempty(ignore_trials)
    n_trials = numel(trl_info.trl_n);
    for f_ix = 1:numel(fields)
        if numel(trl_info.(fields{f_ix}))==n_trials
            trl_info.(fields{f_ix})(ignore_trials) = [];
        end
    end
end
