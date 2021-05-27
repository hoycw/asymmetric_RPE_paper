%% Photodiode Trace Cleaning Parameters: IR87 Run1
% Large drift in the shift over time:
%   plot(evnt(1,:)); Tools --> Basic Fitting --> linear, show equations:
lin_fits   = {[-0.085*evnt.time{1} + 39]};%, ...
%               [-0.05*evnt.time{1} - 43]};
lin_idxs = {[600000:1345000]};
for l_ix = 1:numel(lin_fits)
    evnt.trial{1}(lin_idxs{l_ix}) = evnt.trial{1}(lin_idxs{l_ix})-lin_fits{l_ix}(lin_idxs{l_ix});
end

% Set zero/baseline during a block
bsln_val = 0;

% Record epochs (in sec) with fluctuations that should be set to baseline
% blip between linear fit and bsln_shift chunk; then fixing last blip
bsln_times = {[1344 1345],[1633.6 1633.8]...
    };
% Record epochs (in sec) when baseline has shifted
bsln_shift_times = {[1345 1633.7]};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = [-75];
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
stim_times = {};
stim_yval = [];
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

