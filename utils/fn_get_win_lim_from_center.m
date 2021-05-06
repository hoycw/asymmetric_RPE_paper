function [win_lim] = fn_get_win_lim_from_center(time,win_center,win_len)
%% Find onset and offset of sliding windows given center and width
% INPUTS:
%   time [float] - time vector of the appropriate length and sampling rate
%   win_center [float] - times in seconds to center windows
%   win_len [float] - length in seconds per window
% OUTPUTS:
%   win_lim [int] - N_wins x 2 matrix of indices for [onset offset] of each window

% Check time vector format
if (ndims(time)>2) || ~(size(time,1)==1 || size(time,2)==1)
    error('Too many dimensions in time variable');
end

% Check time vector covers win_center
win_starts = win_center-win_len/2;
win_ends   = win_center+win_len/2;
if min(win_starts)<time(1)-0.00001 || max(win_ends)>time(end)+0.00001
    error('time vector does not cover win_centers!');
end

win_lim = zeros([numel(win_center) 2]);
for win_ix = 1:numel(win_center)
    [~,start_ix] = min(abs(time-win_starts(win_ix)));
    [~,end_ix]   = min(abs(time-win_ends(win_ix)));
    win_lim(win_ix,:) = [start_ix end_ix];
end

missed_ix = setdiff(1:numel(time),min(win_lim(:)):max(win_lim(:)));
if ~isempty(missed_ix)
    fprintf('WARNING!!! %i data points (%3.1f%%) not covered by windows!\n',...
        numel(missed_ix), 100*numel(missed_ix)/numel(time));
end

end