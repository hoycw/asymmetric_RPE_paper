function sorted = fn_sort_labels_alphanum(labels)
%% Sort labels alphanumerically

probes  = cell(size(labels));
lab_num = zeros(size(labels));
for l = 1:numel(labels)
    % Find name of the probe
    probes{l} = labels{l}(regexp(labels{l},'\D'));
    % Find number of each channel
    lab_num(l) = str2num(labels{l}(regexp(labels{l},'\d')));
end
probes_sorted = sort(unique(probes));

% Sort labels by probe then contact number
sorted = {};
for p_ix = 1:numel(probes_sorted)
    cur_lab = labels(strcmp(probes,probes_sorted{p_ix}));
    [~,sort_n_idx] = sort(lab_num(strcmp(probes,probes_sorted{p_ix})));
    sorted = [sorted; cur_lab(sort_n_idx)];
end

end