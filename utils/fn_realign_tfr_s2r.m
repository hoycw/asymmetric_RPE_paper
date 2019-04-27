function [tfr_r] = fn_realign_tfr_s2r(tfr_s,rt,trial_lim_sec)
%% Realign a stimulus-locked TFR to the responses indicated by offsets
% INPUTS:
%   tfr_s [ft struct] - output of ft_freqanalysis to be realigned
%   rt [int array] - n_trialsx1 array of RTs (TIMES in SEC) within that trial
%   trial_lim [int, int] - [start, end] array of TIME (in SEC) indices of epoch around RT
% OUTPUTS:
%   tfr_r [ft_struct] - realigned tfr with correct .powspctrm and .time fields

% Check failure cases
if ~strcmp(tfr_s.dimord,'rpt_chan_freq_time')
    error('Check dimord to be sure trial dimension is first!')
end
if ~isa(tfr_s.time,'double')
    error('tfr_s.time field is not a double. Trials are probably not the same length!');
end

% Compute output powspctrm size
[~, t1_beg_ix] = min(abs(tfr_s.time-(rt(1)+trial_lim_sec(1))));
[~, t1_end_ix] = min(abs(tfr_s.time-(rt(1)+trial_lim_sec(2))));
[~, t1_rt_ix]  = min(abs(tfr_s.time-rt(1)));
new_time = tfr_s.time(t1_beg_ix:t1_end_ix);

% Create new time field
tfr_r = tfr_s;
tfr_r.time = new_time;
tfr_r.time = tfr_r.time-tfr_s.time(t1_rt_ix);
orig_size = size(tfr_s.powspctrm);
tfr_r.powspctrm = NaN([orig_size(1) orig_size(2) orig_size(3) numel(tfr_r.time)]);
for t = 1:size(tfr_s.powspctrm,1)
    [~, beg_ix] = min(abs(tfr_s.time-(rt(t)+trial_lim_sec(1))));
    [~, end_ix] = min(abs(tfr_s.time-(rt(t)+trial_lim_sec(2))));
    if numel(beg_ix:end_ix)~=numel(tfr_r.time)
        warning('sampling rate causes slight misalignments in time!');
        end_ix = end_ix-1;
    end
    tfr_r.powspctrm(t,:,:,:) = tfr_s.powspctrm(t,:,:,beg_ix:end_ix);
end

end