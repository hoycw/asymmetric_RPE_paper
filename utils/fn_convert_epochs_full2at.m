function cut_epochs = fn_convert_epochs_full2at(epochs,analysis_time,full_filename,keep_partial)
%% Convert epochs from original sample indices to trimmed analysis_time indices by matching data.time indices
% INPUTS:
%   epochs [int array] - sample indices in original time series (not yet cut)
%   analysis_time [cell array of tuples] - contains [start stop] times in SECONDS of valid segments of data
%   full_filename [str] - full path to the original time series data
%       probably "SBJ_preclean.mat"
%   keep_partial [0/1] - flag to handle cases where epochs stradle segment borders
%       0: remove epochs without FULL overlap with valid segments
%       1: trim epoch edges to include only valid points
%       NOTE: completely non-overlapping epochs will be removed in both cases!

addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Load full file
load(full_filename);

%% Convert epoch samples into times
orig_times = NaN(size(epochs));
for ep_ix = 1:size(epochs,1)
    orig_times(ep_ix,:) = [data.time{1}(epochs(ep_ix,1)) data.time{1}(epochs(ep_ix,2))];
end

%% Cut full data to analysis_time
for seg_ix = 1:length(analysis_time)
    epoch_len{seg_ix} = diff(analysis_time{seg_ix});
    cfg_trim = [];
    cfg_trim.latency = analysis_time{seg_ix};
    
    data_pieces{seg_ix} = ft_selectdata(cfg_trim, data);
end
% Stitch them back together
new_time = data_pieces{1}.time{1};  % this keeps the old times without making them continuous!
if length(analysis_time)>1
    for seg_ix = 2:length(analysis_time)
        new_time = horzcat(new_time,data_pieces{seg_ix}.time{1});
    end
end

%% Find new indices of original times
cut_epochs = NaN(size(epochs));
for ep_ix = 1:size(epochs,1)
    new_ix1 = find(new_time==orig_times(ep_ix,1));
    new_ix2 = find(new_time==orig_times(ep_ix,2));
    % If overlapping at all with valid segment, adjust epoch (otherwise, leave as NaN!)
    if ~isempty(new_ix1) || ~isempty(new_ix2)
        if ~isempty(new_ix1) && ~isempty(new_ix2)   % Both are good, adjust them
            cut_epochs(ep_ix,:) = [new_ix1, new_ix2];
        elseif keep_partial % Only one is good AND keep_partial==1, adjust edge to fit in valid times
            if isempty(new_ix1)
                cut_epochs(ep_ix,1) = find(new_time>orig_times(ep_ix,1),1,'first');
                cut_epochs(ep_ix,2) = new_ix2;
            else
                assert(isempty(new_ix2));
                cut_epochs(ep_ix,1) = new_ix1;
                cut_epochs(ep_ix,2) = find(new_time<orig_times(ep_ix,2),1,'last');
            end
        end
%     else
%         % If before segment starts, do nothing
%         if orig_times(ep_ix,2)<new_time(1)
%             continue;
%         elseif epochs(ep_ix,1)>seg_end
%             seg_ix    = seg_ix+1;
%             % If completely past segment, advance segment or end loop
%             if seg_ix>numel(analysis_samples)
%                 break;
%             else
%                 seg_start = analysis_samples{seg_ix}(1);
%                 seg_end   = analysis_samples{seg_ix}(2);
%                 offset    = analysis_seg_len{seg_ix-1};
%             end
%         end
    end
end

%% Get rid of NaN (non-matching) epochs
cut_epochs(isnan(cut_epochs(:,1)),:) = [];

%% Report results
if size(cut_epochs,1)==size(epochs,1)
    fprintf('SUCCESS! All %i epochs were adjusted.\n',size(epochs,1));
else
    fprintf('WARNING!!! Started with %i epochs, only %i were adjusted!\n',...
        size(epochs,1),size(cut_epochs,1));
end

end