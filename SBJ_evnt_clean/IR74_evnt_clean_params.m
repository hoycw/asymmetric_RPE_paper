%% Photodiode Trace Cleaning Parameters: IR74
% Set zero/baseline during a block
bsln_val = 35310;
evnt.trial{1}(end) = bsln_val;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0 2],... % feedback on fullvis trial 2
    [252 576],... % delete B1T21+22 through B1T28 with 2 pauses, resume at B1T29
    [1538 1539]... % blip at end
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

