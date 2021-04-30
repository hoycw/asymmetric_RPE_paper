%% Photodiode Trace Cleaning Parameters: IR84
% Mark trials to ignore e.g., interruptions
%   first 5 full vis must've been cut off in clinical extraction
ignore_trials = [1:5];

% Large drift in the shift over time:
%   plot(evnt(1,:)); Tools --> Basic Fitting --> linear, show equations:
lin_fits   = {[-0.12*evnt.time{1} + 45], ...
              [-0.05*evnt.time{1} - 43]};
lin_idxs = {[67850:450000], [450001:720000]};
for l_ix = 1:numel(lin_fits)
    evnt.trial{1}(lin_idxs{l_ix}) = evnt.trial{1}(lin_idxs{l_ix})-lin_fits{l_ix}(lin_idxs{l_ix});
end

% Set zero/baseline during a block
bsln_val = 0;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
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

