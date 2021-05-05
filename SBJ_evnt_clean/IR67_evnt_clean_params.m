%% Photodiode Trace Cleaning Parameters: IR67
% Set zero/baseline during a block
bsln_val = -318;

% fix initial baseline (different value than ending) 
[~,time_ix] = min(abs(evnt.time{1}-48));
evnt.trial{1}(1:time_ix) = -218;

% get rid of last value
evnt.trial{1}(end) = bsln_val;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {[1255 1380]...
    };
% Record epochs (in sec) when baseline has shifted
bsln_shift_times = {};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = [];
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
stim_times = {};
stim_yval = [];
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

