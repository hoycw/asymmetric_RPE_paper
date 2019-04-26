function [bhv] = fn_load_behav_csv(csv_fname,ignore_trials)
%% Load trial info csv exported from python
%   Also converts from python 0-based to MATLAB 1-based indexing
% INPUTS:
%   csv_fname [str] - full file path to behavioral csv
%   ignore_trials [int array] - list of trial numbers to be removed (MATLAB indexing)
% OUTPUTS:
%   bhv [struct] - contains many fields of info on the trials

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


% Open file
fprintf('\tReading behavioral csv file: %s\n',csv_fname);
bhv_file = fopen(csv_fname);

% Read field names
%   Fields: Total_Trial, Block, Condition, Hit, RT, Timestamp,...
%           Tolerance, Trial, Score, ITI, ITI type
% fields = {'trl_n','blk','cond','hit','rt','time',...
%             'tol','blk_trl','score','ITI','ITI_type'};
bhv_fields = textscan(bhv_file, '%s %s %s %s %s %s %s %s %s %s %s', 1, 'Delimiter', ',');

% Read data
%   orig version: formatspec = '%d%d%s%d%f%f%f%d%d%f%s';
bhv_data = textscan(bhv_file,'%d %d %s %d %f %f %f %d %d %f %f',...
                'Delimiter',',','HeaderLines',1);
fclose(bhv_file);

% Get list of good trials
n_trials = size(bhv_data{1},1);
fprintf('\t\tFound %d trials in log file\n', n_trials);
if ~isempty(ignore_trials)
    fprintf('\t\tIgnoring %d trials\n', numel(ignore_trials));
    good_trials = setdiff(1:n_trials,ignore_trials);
else
    good_trials = 1:n_trials;
end

% Add data to struct, convert field names with spaces to underscore
for f_ix = 1:numel(bhv_fields)
    bhv_fields{f_ix} = strrep(bhv_fields{f_ix}{1},' ','_');
    if numel(bhv_data{f_ix})==n_trials
        bhv.(bhv_fields{f_ix}) = bhv_data{f_ix}(good_trials);
    else
        bhv.(bhv_fields{f_ix}) = bhv_data{f_ix};
    end
end

% Convert from python to MATLAB indexing
bhv.Total_Trial = bhv.Total_Trial+1;
bhv.Block       = bhv.Block+1;
bhv.Trial       = bhv.Trial+1;

end
