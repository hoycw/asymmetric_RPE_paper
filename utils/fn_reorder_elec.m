function [ordered_elec] = fn_reorder_elec(elec, labels)
%% Re-order an elec struct to match order of labels provided
% If labels is empty, sort alphabetically
if isempty(labels)
    labels = fn_sort_labels_alphanum(elec.label);
end

ordered_elec = elec;
% Order them to match labels provided
[matches,order_idx] = ismember(labels,elec.label);
if sum(matches)~=numel(labels)
    fprintf('Mismatched electrodes not found in elec:\n');
    disp(labels(matches==0));
    error('Mismatch between labels provided and elec.label!');
else
    fields = fieldnames(elec);
    for f = 1:numel(fields)
        if size(eval(['elec.' fields{f}]),1) == numel(elec.label)
            eval(['ordered_elec.' fields{f} ' = elec.' fields{f} '(order_idx(matches==1),:);']);
        elseif size(eval(['elec.' fields{f}]),1) > numel(elec.label)
            warning(['elec field "' fields{f} '" has more elements than elec.label, and will not be re-ordered!!!']);
        end 
    end 
end

end
