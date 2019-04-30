%% Photodiode Trace Cleaning Parameters: CP24 
% Mark trials to ignore e.g., interruptions
ignore_trials = [];

% Fix drift over first half of recording
x = [1:find(evnt.time{1}==867)];
y = -0.000042201*x + 283.93;
evnt.trial{1}(x) = evnt.trial{1}(x)-y;
% Remove offset in remainder
evnt.trial{1}(x(end)+1:end) = evnt.trial{1}(x(end)+1:end)-249;

% Set zero/baseline during a block
bsln_val = 0;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0.0 157.0]...%restart
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

