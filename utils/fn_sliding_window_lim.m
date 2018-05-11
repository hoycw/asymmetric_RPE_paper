function [win_lim] = fn_sliding_window_lim(time,win_len,win_step)
%% Find onset and offset of sliding windows
% INPUTS:
%   trial [float] - one sample time series of the appropriate length
%   win_len [int] - number of data points per window
%   win_step [int] - number of data points to advance each window
% OUTPUTS:
%   win_lim [int] - N_wins x 2 matrix of indices for [onset offset] of each window

if (ndims(time)>2) || ~(size(time,1)==1 || size(time,2)==1)
    error('Too many dimensions in time variable');
end

win_lim = zeros([ceil((numel(time)-win_len-1)/win_step) 2]);
win_lim(1,:) = [1 win_len];
for win_ix = 2:size(win_lim,1)
    win_lim(win_ix,:) = [win_lim(win_ix-1,1)+win_step win_lim(win_ix-1,2)+win_step];
end

if win_lim(end,2)<numel(time)
    fprintf('WARNING!!! %i data points (%3.1f%%) not covered by windows!\n',...
        numel(time(win_lim(end,2)+1:end)),...
        100*numel(time(win_lim(end,2)+1:end))/numel(time));
end

end