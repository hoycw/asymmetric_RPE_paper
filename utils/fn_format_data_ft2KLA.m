function [data_ecog, header_ecog] = fn_format_data_ft2KLA(data_ft)
%% Save ft to data_ecog format

data_ecog = data_ft.trial{:};

header_ecog.sample_rate = data_ft.fsample;
header_ecog.orig_n_channels = size(data_ft.label,1);
header_ecog.n_samples = size(data_ecog,2);
header_ecog.length_in_seconds = size(data_ft.trial{:},2)/data_ft.fsample;
header_ecog.original_sample_rate = data_ft.fsample;
header_ecog.channel_labels = {data_ft.label{:}};

header_ecog.n_channels = size(data_ecog,1);
% header_ecog.orig_channel_n = analysis_channels;
% header_ecog.raw_file = raw_filename;

end