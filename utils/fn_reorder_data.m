function [ordered_data] = fn_reorder_data(data, labels)
%% Re-order an data struct to match order of labels provided
% If labels is empty, sort alphabetically
if isempty(labels)
    labels = fn_sort_labels_alphanum(data.label);
end

ordered_data = data;
% Order them to match labels provided
[matches,order_idx] = ismember(labels,data.label);
if sum(matches)~=numel(labels)
    fprintf('Mismatched electrodes from input labels not found in data:\n');
    disp(label(matches==0));
    error('Mismatch between labels provided and data.label!');
else
    fields = fieldnames(data);
    for f = 1:numel(fields)
        if size(eval(['data.' fields{f}]),1) == numel(data.label)
            eval(['ordered_data.' fields{f} ' = data.' fields{f} '(order_idx(matches==1),:);']);
        elseif size(eval(['data.' fields{f}]),1) > numel(data.label)
            warning(['data field "' fields{f} '" has more elements than data.label, and will not be re-ordered!!!']);
        end 
    end 
    for t = 1:numel(data.trial)
        ordered_data.trial{t} = data.trial{t}(order_idx(matches==1),:);
    end
end

end
