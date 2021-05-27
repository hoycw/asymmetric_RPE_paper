%% Photodiode Trace Cleaning Parameters: IR94 Run2
% Large drift in the shift over time:
%   plot(evnt(1,:)); Tools --> Basic Fitting --> linear, show equations:
lin_fits   = {[-0.06*evnt.time{1} -5]};%, ...
lin_idxs = {[700000:1385001]};
for l_ix = 1:numel(lin_fits)
    % -50 to account for baseline
    evnt.trial{1}(lin_idxs{l_ix}) = evnt.trial{1}(lin_idxs{l_ix})-lin_fits{l_ix}(lin_idxs{l_ix})-50;
end

% Set zero/baseline during a block
bsln_val = -35;
% Fix initial offset
evnt.trial{1}(1:20000) = bsln_val;
evnt.trial{1}(end) = bsln_val;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {[640 710],[1269.5 1385]...
    };
% Record epochs (in sec) when baseline has shifted
bsln_shift_times = {[710 1385]};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = [-20];
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
stim_times = {};
stim_yval = [];
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

