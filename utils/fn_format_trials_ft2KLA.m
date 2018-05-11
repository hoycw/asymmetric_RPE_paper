function [trials, header_ecog] = fn_format_trials_ft2KLA(trials_ft)
%% Save ft to data_ecog format

trials = NaN([size(trials_ft.trial{1},1) numel(trials_ft.trial) size(trials_ft.trial{1},2)]);
for t = 1:numel(trials_ft.trial)
    trials(:,t,:) = trials_ft.trial{t};
end

header_ecog.sample_rate = trials_ft.fsample;
header_ecog.n_samples = size(trials,1)*size(trials,2)*size(trials,3);
header_ecog.channel_labels = {trials_ft.label{:}};
header_ecog.n_channels = size(trials,1);

end