function [data_clean] = fn_cleanline(data,line_noise_freqs)
%% Run Tim Mullen's cleanline on this data
% Some important comments from Readme.txt:
%   -If cleaning continuous, un-epoched data, then you may wish to use sliding
%       windows of 3-4 seconds with 50% overlap.

addpath(genpath('/home/knight/hoycw/PRJ_Error/scripts/utils/'));

cleaned = cleanline2(data.trial{1}, data.fsample,...
    'LineFrequencies',  line_noise_freqs, ...   % default is 60, 120
    'ScanForLines',     1, ...                  % finds exact line freq around given value
    'LineAlpha',        0.01, ...               % default = 0.01
    'Bandwidth',        3, ...                  % default = 1
    'ChanCompIndices',  1:numel(data.label), ...
    'SlidingWinLength', 4.0, ...                % in sec, default = 4
    'SlidingWinStep',   2.0, ...                % in sec, default = 4 (no overlap)
    'SmoothingFactor',  100, ...                % default=100; 1=linear, Inf=no smooth
    'PaddingFactor',    2, ...                  % default = 2; KLA used 1
    'ComputeSpectralPower',0,...                % NOT POSSIBLE in cleanline2
    'PlotFigures', 1 ...                        % NOT POSSIBLE in cleanline2- maybe alternative to CompSpecPow
    );

% Update header information
data_clean = data;
data_clean.trial{1} = cleaned;
data_clean.cfg.cleanlinefilter = 'yes';
data_clean.cfg.notch_freq = line_noise_freqs;

end
