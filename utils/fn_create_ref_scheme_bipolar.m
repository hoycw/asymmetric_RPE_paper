function [new_labels, weights, left_out_ch] = fn_create_ref_scheme_bipolar(labels)
%% Create a re-reference scheme with bipolar pairs from a list of labels
%   Goes from smallest to largest (ELEC1-ELEC2, ELEC2-ELEC3, etc.)
%   Pairs made only from contiguous numbers
%
% new_labels = cell array of strings with combined labels
%   can be used with ft_preprocessing as cfg_rereference.montage.labelnew
% weights = [length(new_labels),length(labels)] matrix of weights to combine elecs
%   can be used with ft_preprocessing as cfg_rereference.montage.tra
% left_out_ch = labels of channels that aren't included (don't have contiguous pair)

% Find name of the probe
probe = labels{1}(regexp(labels{1},'\D'));

% Find number of each channel
lab_num = zeros(size(labels));
for lab_ix = 1:length(labels)
    if ~strcmp(probe,labels{lab_ix}(regexp(labels{lab_ix},'\D')))
        error('ERROR!!! Mixed probe names in labels!');
    end
    lab_num(lab_ix) = str2num(labels{lab_ix}(regexp(labels{lab_ix},'\d')));
end

% Combine pairs into labels and weights
weights = zeros(length(find(diff(lab_num)==1)),numel(labels));
new_labels = {};
weight_row = 0;
for lab_ix = 2:numel(labels)
    if lab_num(lab_ix)-lab_num(lab_ix-1)==1
        new_labels = {new_labels{:} strcat(labels{lab_ix-1},'-',num2str(lab_num(lab_ix)))};
        weight_row = weight_row+1;
        weights(weight_row,lab_ix-1) = 1;
        weights(weight_row,lab_ix)   = -1;
    end
end

% Find any channels left out
left_out_ch = {};
for col_ix = 1:size(weights,2)
    if length(find(weights(:,col_ix)==0))==size(weights,1)
        left_out_ch = {left_out_ch{:} labels{col_ix}};
    end
end

% Warn user of left out channels
if ~isempty(left_out_ch)
    fprintf('===========================================================\n');
    fprintf('!!!WARNING!!! %i channels without pairs!\nLeft out channels: ',numel(left_out_ch));
    for ix = 1:numel(left_out_ch)
        fprintf('%s',left_out_ch{ix});
    end
    fprintf('\n');
    fprintf('===========================================================\n');
end

end
