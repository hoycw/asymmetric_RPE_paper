%% Photodiode Trace Cleaning Parameters: IR94 Run2
% Large drift in the shift over time:
%   plot(evnt(1,:)); Tools --> Basic Fitting --> linear, show equations:

% Fix initial offset
evnt.trial{1}(1:20000) = -100;

%% Set zero/baseline during a block
bsln_val = -200;
evnt_val = 440;

%% Record conditional within trial corrections
cond_times  = {...
    {[440 479.5],evnt_val,-15},...      % B1 T1-11 SR events
    {[479.5 481],evnt_val,-200},...     % B1 T11 F events
    {[481 489.8],evnt_val,-50},...   % B1 T12 to T14 onset
    {[487.8 488.6],evnt_val,-200},...   % B1 T13 FB (drift)
    {[491.3 506.7],evnt_val,-200},...   % B1 T14 FB + T15:T18 events
    {[1657 1659],evnt_val,-100},...     % B4 T75 FB events
    };
%     {[440 479.5],bsln_val,-15},...      % B1 T1-11 SR baselines
%     {[491.2 506.7],bsln_val,-200},...   % B1 T14 R + T15:T18 baseline
%     {[1657 1700],bsln_val,100}...       % end of recording baseline

for cond_ix = 1:numel(cond_times)
    if numel(cond_times{cond_ix})~=3 || numel(cond_times{cond_ix}{1})~=2
        error(['check cond_times ix = ' num2str(cond_ix)]);
    end
    epoch_ix = floor(cond_times{cond_ix}{1}(1)*evnt.fsample):floor(cond_times{cond_ix}{1}(2)*evnt.fsample);
    epoch_ix(epoch_ix<1) = [];
    cond_idx = evnt.trial{1}(epoch_ix)>cond_times{cond_ix}{3};
    evnt.trial{1}(epoch_ix(cond_idx)) = evnt_val;
    
%     if cond_set(cond_ix)==evnt_val
%         % greater than
%     elseif cond_set(cond_ix)==bsln_val
%         % less than
%     else
%         error('why set to something besides event or baseline levels?');
%     end
end

%% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0 440],...     % initial mess (full vis, both trainings)
    [497 498],...   % B1 T16 response period drifted up
    [1680 1700]...  % end of NLX, set to 0
    };
% zero out drifts
for zero_ix = 1:length(bsln_times)
    epoch_ix = floor(bsln_times{zero_ix}(1)*evnt.fsample):floor(bsln_times{zero_ix}(2)*evnt.fsample);
    epoch_ix(epoch_ix<1) = [];
    evnt.trial{1}(epoch_ix) = bsln_val;
end

% Fix blip at end
evnt.trial{1}(end) = bsln_val;

%% Record epochs (in sec) when baseline has shifted
bsln_shift_times = {};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = [];
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end
% Correct baseline shift
for shift_ix = 1:length(bsln_shift_times)
    epoch_ix = floor(bsln_shift_times{shift_ix}(1)*evnt.fsample):floor(bsln_shift_times{shift_ix}(2)*evnt.fsample);
    epoch_ix(epoch_ix<1) = [];
    evnt.trial{1}(epoch_ix) = evnt.trial{1}(epoch_ix) - bsln_shift_val(shift_ix);
end

%% Record within trial corrections
stim_times = {[484 484.7]};
stim_yval = [evnt_val];
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

% level out stimulus periods
for stim_ix = 1:length(stim_times)
    epoch_ix = floor(stim_times{stim_ix}(1)*evnt.fsample):floor(stim_times{stim_ix}(2)*evnt.fsample);
    epoch_ix(epoch_ix<1) = [];
    evnt.trial{1}(epoch_ix) = stim_yval(stim_ix);
end

%% Final binarization
binary_thresh = 0;
evnt.trial{1}(evnt.trial{1}<binary_thresh)  = bsln_val;
evnt.trial{1}(evnt.trial{1}>=binary_thresh) = evnt_val;

