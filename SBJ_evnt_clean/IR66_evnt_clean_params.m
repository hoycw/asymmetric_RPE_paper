%% Photodiode Trace Cleaning Parameters: IR68
% Set zero/baseline during a block
bsln_val = 21550;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0.0 10.0],...%blip at start
    [750.0 760.0],...% blip between blocks
    [1338.0 1346.09]...
    };
evnt.trial{1}(end) = bsln_val;

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

